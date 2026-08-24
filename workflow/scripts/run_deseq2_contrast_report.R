#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(EnhancedVolcano)
  library(rmarkdown)
})

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

default_if_null <- function(value, default) {
  if (is.null(value)) {
    return(default)
  }
  value
}

as_bool <- function(value, name) {
  if (is.logical(value) && length(value) == 1 && !is.na(value)) {
    return(value)
  }
  if (is.numeric(value) && length(value) == 1 && !is.na(value)) {
    return(value != 0)
  }
  if (is.character(value) && length(value) == 1) {
    lowered <- tolower(trimws(value))
    if (lowered %in% c("true", "1", "yes")) {
      return(TRUE)
    }
    if (lowered %in% c("false", "0", "no")) {
      return(FALSE)
    }
  }
  stopf("%s must be a boolean value. Found: %s", name, paste(value, collapse = ","))
}

# Names of every synthetic/computed column that gets attached to a count
# matrix's metadata columns downstream -- by run_deseq2_matrix() below via
# dplyr::mutate() into result_tbl, and by deseq2_contrast_report.Rmd's
# make_volcano()/make_ma() via a further mutate() on that same table. A
# metadata column sharing one of these names would be silently overwritten
# by mutate() rather than erroring, so this list backs an explicit collision
# check instead.
RESERVED_RESULT_COLUMNS <- c(
  "comparison", "group1", "group2",
  "baseMean", "log2FoldChange_raw", "log2FoldChange_shrunk", "shrinkage_applied",
  "lfcSE", "stat", "pvalue", "padj", "adjusted_pvalue", "mean_group1", "mean_group2",
  "significant", "direction", "neg_log10_padj", "feature_label"
)

split_delimited <- function(value, delimiter = ";;;") {
  value <- trimws(default_if_null(value, ""))
  if (identical(value, "")) {
    return(character(0))
  }
  pieces <- strsplit(value, delimiter, fixed = TRUE)[[1]]
  trimws(pieces[nzchar(trimws(pieces))])
}

# Human-readable names for the fixed set of repeat classes defined by
# REPEAT_TYPES in workflow/rules/deseq2.smk (matched on their lowercased key).
REPEAT_CLASS_LABELS <- c(
  sine_alu = "SINE Alu",
  sine_mir = "SINE MIR",
  line_l1 = "LINE L1",
  line_l2 = "LINE L2",
  ltr = "LTR",
  other_repeat_elements = "Other Repeat Elements"
)

# Maps an internal DESEQ2_AVAILABLE_MATRICES key (e.g. "trna", "pol3_t1",
# "repeat_sine_alu", "viral_KT899744.1") to a human-readable report section
# title. Only "viral_*" keys are actual virus genomes -- everything else is
# host-associated (tRNA/Pol III/repeat-element/rRNA), despite historically
# being plumbed through CLI flags named --virus-*.
category_display_name <- function(feature_type) {
  if (identical(feature_type, "host")) {
    return("Host")
  }
  if (identical(feature_type, "trna")) {
    return("tRNA")
  }
  if (identical(feature_type, "rrna")) {
    return("rRNA")
  }
  if (grepl("^pol3_t[0-9]+$", feature_type)) {
    return(sprintf("Pol III (%s)", toupper(sub("^pol3_", "", feature_type))))
  }
  if (grepl("^viral_", feature_type)) {
    return(sprintf("Virus: %s", sub("^viral_", "", feature_type)))
  }
  if (grepl("^repeat_", feature_type)) {
    repeat_key <- sub("^repeat_", "", feature_type)
    repeat_name <- REPEAT_CLASS_LABELS[[repeat_key]]
    if (is.null(repeat_name)) {
      repeat_name <- tools::toTitleCase(gsub("_", " ", repeat_key))
    }
    return(sprintf("Repeat Element: %s", repeat_name))
  }
  tools::toTitleCase(gsub("_", " ", feature_type))
}

pick_feature_label_column <- function(df) {
  candidates <- c("gene_name", "gene_id", "bin_id", "chrom")
  for (candidate in candidates) {
    if (candidate %in% colnames(df)) {
      values <- as.character(df[[candidate]])
      if (any(nzchar(values))) {
        return(candidate)
      }
    }
  }
  NA_character_
}

safe_feature_labels <- function(df, preferred_column) {
  if (!is.na(preferred_column) && preferred_column %in% colnames(df)) {
    labels <- as.character(df[[preferred_column]])
    labels[!nzchar(labels)] <- NA_character_
  } else {
    labels <- rep(NA_character_, nrow(df))
  }

  if ("bin_id" %in% colnames(df)) {
    labels[is.na(labels)] <- as.character(df$bin_id[is.na(labels)])
  }

  has_coords <- all(c("chrom", "start", "end") %in% colnames(df))
  if (has_coords) {
    coord_labels <- paste0(df$chrom, ":", df$start, "-", df$end)
    labels[is.na(labels)] <- coord_labels[is.na(labels)]
  }

  labels[is.na(labels)] <- paste0("feature_", seq_len(sum(is.na(labels))))
  labels
}

save_enhanced_volcano_png <- function(result_tbl, feature_labels, png_path, comparison,
                                      group1, group2, alpha, lfc_threshold, label_top_n) {
  padj_safe <- ifelse(is.na(result_tbl$adjusted_pvalue), 1, result_tbl$adjusted_pvalue)

  sig_mask <- !is.na(result_tbl$adjusted_pvalue) &
    result_tbl$adjusted_pvalue < alpha &
    !is.na(result_tbl$log2FoldChange_shrunk) &
    abs(result_tbl$log2FoldChange_shrunk) >= lfc_threshold

  select_lab <- rep("", nrow(result_tbl))
  if (sum(sig_mask) > 0) {
    sig_idx <- which(sig_mask)
    ord <- order(result_tbl$adjusted_pvalue[sig_idx], -abs(result_tbl$log2FoldChange_shrunk[sig_idx]))
    top_idx <- sig_idx[ord][seq_len(min(label_top_n, length(sig_idx)))]
    select_lab[top_idx] <- feature_labels[top_idx]
  }

  p <- EnhancedVolcano::EnhancedVolcano(
    toptable        = data.frame(log2FC = result_tbl$log2FoldChange_shrunk,
                                 padj   = padj_safe,
                                 stringsAsFactors = FALSE),
    lab             = select_lab,
    x               = "log2FC",
    y               = "padj",
    title           = sprintf("%s  |  %s vs %s", comparison, group1, group2),
    subtitle        = bquote(alpha == .(alpha) ~ "| LFC ==" ~ .(lfc_threshold)),
    pCutoff         = alpha,
    FCcutoff        = lfc_threshold,
    xlab            = bquote(Shrunken ~ log[2] ~ fold ~ change),
    ylab            = bquote(-log[10] ~ p[adj]),
    legendPosition  = "right",
    drawConnectors  = TRUE,
    widthConnectors = 0.5,
    colConnectors   = "grey30",
    pointSize       = 1.5,
    labSize         = 3.0,
    max.overlaps    = Inf
  )

  ggplot2::ggsave(filename = png_path, plot = p,
                  width = 10, height = 8, dpi = 150, units = "in")
  invisible(png_path)
}

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg) || !is.list(cfg)) {
    stopf("Config file did not parse into a mapping: %s", path)
  }
  cfg
}

# deseq2.skip_features is matched against the literal, lowercased internal
# category key for each matrix (e.g. "trna", "pol3_t1", "repeat_sine_alu",
# "viral_<accession>") -- it is not a grouping/alias mechanism, so entries
# like "virus" or "ALU" never match anything and silently skip nothing. Warn
# loudly (rather than abort) since an unmatched entry is inert, not
# destructive: the run still produces correct output, just not the reduced
# scope the user intended.
warn_unmatched_skip_features <- function(skip_features, available_keys) {
  if (length(skip_features) == 0) {
    return(invisible(NULL))
  }
  available_lower <- tolower(available_keys)
  unmatched <- skip_features[!(skip_features %in% available_lower)]
  if (length(unmatched) > 0) {
    cat(sprintf(
      paste(
        "\nWARNING: deseq2.skip_features contains %s that match no available",
        "category for this run and will skip nothing: %s",
        "Available categories: %s",
        "skip_features matches literal category keys, not groupings like",
        "'virus' or 'ALU' -- e.g. use 'viral_<accession>' for a specific",
        "virus, or list each repeat/Pol III class you want skipped.\n",
        sep = "\n"
      ),
      if (length(unmatched) == 1) "an entry" else "entries",
      paste(unmatched, collapse = ", "),
      paste(sort(available_lower), collapse = ", ")
    ))
  }
  invisible(NULL)
}

validate_design_factors <- function(design_factors) {
  if (is.null(design_factors)) {
    stopf("deseq2.design_factors is not set.")
  }
  if (is.character(design_factors) && length(design_factors) == 1) {
    design_factors <- trimws(design_factors)
  }
  if (identical(design_factors, "group")) {
    return(invisible(NULL))
  }
  if (is.list(design_factors)) {
    design_factors <- unlist(design_factors, use.names = FALSE)
  }
  if (!identical(as.character(design_factors), "group")) {
    stopf("Only deseq2.design_factors = ['group'] is implemented in this workflow.")
  }
}

build_deseq2_config <- function(cfg) {
  dcfg <- default_if_null(cfg$deseq2, list())
  list(
    alpha = as.numeric(default_if_null(dcfg$alpha, 0.05)),
    lfc_threshold = as.numeric(default_if_null(dcfg$lfc_threshold, 1.0)),
    min_total_count = as.integer(default_if_null(dcfg$min_total_count, 1)),
    size_factor_type = as.character(default_if_null(dcfg$size_factor_type, "poscounts")),
    fit_type = as.character(default_if_null(dcfg$fit_type, "parametric")),
    shrink_type = as.character(default_if_null(dcfg$shrink_type, "apeglm")),
    skip_features = tolower(unlist(default_if_null(dcfg$skip_features, list()), use.names = FALSE)),
    min_features_per_matrix = as.integer(default_if_null(dcfg$min_features_per_matrix, 50)),
    p_adjust_method = as.character(default_if_null(dcfg$p_adjust_method, "BH")),
    cooks_cutoff = as_bool(default_if_null(dcfg$cooks_cutoff, TRUE), "deseq2.cooks_cutoff"),
    independent_filtering = as_bool(
      default_if_null(dcfg$independent_filtering, TRUE),
      "deseq2.independent_filtering"
    ),
    vst_blind = as_bool(default_if_null(dcfg$vst_blind, TRUE), "deseq2.vst_blind"),
    report_top_n = as.integer(default_if_null(dcfg$report_top_n, 50)),
    report_label_top_n = as.integer(default_if_null(dcfg$report_label_top_n, 15)),
    report_max_table_rows = as.integer(default_if_null(dcfg$report_max_table_rows, 50000)),
    html_self_contained = as_bool(
      default_if_null(dcfg$html_self_contained, TRUE),
      "deseq2.html_self_contained"
    )
  )
}

read_manifest <- function(path, group1, group2) {
  manifest <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  required_cols <- c("sampleName", "groupName", "batch")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) {
    stopf("Manifest is missing required columns: %s", paste(missing_cols, collapse = ", "))
  }

  selected <- manifest %>%
    dplyr::filter(.data$groupName %in% c(group1, group2)) %>%
    dplyr::mutate(
      sampleName = as.character(.data$sampleName),
      groupName = as.character(.data$groupName),
      batch = as.character(.data$batch)
    )

  if (nrow(selected) == 0) {
    stopf("No manifest samples matched groups '%s' and '%s'.", group1, group2)
  }

  if (!all(c(group1, group2) %in% selected$groupName)) {
    stopf("Manifest selection did not retain both groups '%s' and '%s'.", group1, group2)
  }

  list(
    manifest = manifest,
    selected = selected
  )
}

prepare_counts <- function(matrix_df, all_manifest_samples, selected_samples, min_total_count, feature_type) {
  metadata_cols <- setdiff(colnames(matrix_df), all_manifest_samples)
  reserved_collisions <- intersect(metadata_cols, RESERVED_RESULT_COLUMNS)
  if (length(reserved_collisions) > 0) {
    stopf(
      "Matrix metadata column(s) collide with reserved DESeq2 result column name(s) for %s: %s. Rename the matrix column(s) to avoid silent data loss.",
      feature_type, paste(reserved_collisions, collapse = ", ")
    )
  }
  missing_samples <- setdiff(selected_samples, colnames(matrix_df))
  if (length(missing_samples) > 0) {
    stopf(
      "Matrix is missing required sample columns: %s",
      paste(missing_samples, collapse = ", ")
    )
  }

  counts_df <- matrix_df[, selected_samples, drop = FALSE]
  counts_mat <- as.matrix(counts_df)
  storage.mode(counts_mat) <- "double"
  if (any(is.na(counts_mat))) {
    stopf("Matrix contains NA count values for selected samples.")
  }
  if (any(counts_mat < 0)) {
    stopf("Matrix contains negative counts, which DESeq2 does not support.")
  }
  if (any(abs(counts_mat - round(counts_mat)) > 1e-8)) {
    stopf("Matrix contains non-integer counts, which DESeq2 does not support.")
  }
  counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"

  keep <- rowSums(counts_mat) >= min_total_count
  if (!any(keep)) {
    stopf(
      "No features passed min_total_count=%d for the requested comparison.",
      min_total_count
    )
  }

  metadata_df <- as.data.frame(matrix_df[keep, metadata_cols, drop = FALSE], stringsAsFactors = FALSE)
  row_ids <- NULL
  if ("bin_id" %in% colnames(metadata_df)) {
    row_ids <- as.character(metadata_df$bin_id)
  } else if ("gene_id" %in% colnames(metadata_df)) {
    row_ids <- as.character(metadata_df$gene_id)
  }
  if (is.null(row_ids)) {
    row_ids <- paste0(feature_type, "_", seq_len(nrow(metadata_df)))
  }
  row_ids[is.na(row_ids) | !nzchar(row_ids)] <- paste0(feature_type, "_", seq_len(sum(is.na(row_ids) | !nzchar(row_ids))))
  row_ids <- make.unique(row_ids)

  rownames(metadata_df) <- row_ids
  counts_mat <- counts_mat[keep, , drop = FALSE]
  rownames(counts_mat) <- row_ids

  list(
    metadata_df = metadata_df,
    counts_mat = counts_mat,
    n_features = nrow(metadata_df)
  )
}

run_deseq2_matrix <- function(
  matrix_path,
  output_path,
  manifest,
  selected_manifest,
  comparison,
  group1,
  group2,
  feature_type,
  cfg
) {
  feature_type_lower <- tolower(feature_type)

  if (feature_type_lower %in% cfg$skip_features) {
    skip_reason <- sprintf(
      "Feature type '%s' is in the skip list (configured: %s)",
      feature_type, paste(cfg$skip_features, collapse = ", ")
    )
    cat(sprintf("\n⏭️  ANALYSIS SKIPPED: %s - %s\n%s\n", comparison, feature_type, skip_reason))
    readr::write_lines(
      sprintf(
        "# ANALYSIS SKIPPED: %s - %s\n# Reason: Feature type in skip list\n# Skip list: %s",
        comparison, feature_type, paste(cfg$skip_features, collapse = ", ")
      ),
      output_path
    )
    return(list(
      skipped = TRUE, skip_type = "skip_list", reason = skip_reason,
      feature_type = feature_type, display_name = category_display_name(feature_type),
      n_features = NA_integer_
    ))
  }

  matrix_df <- suppressMessages(readr::read_tsv(matrix_path, show_col_types = FALSE, progress = FALSE))
  if (nrow(matrix_df) == 0) {
    stopf("Matrix has no rows: %s", matrix_path)
  }

  all_manifest_samples <- as.character(manifest$sampleName)
  selected_samples <- as.character(selected_manifest$sampleName)
  prepared <- prepare_counts(
    matrix_df = matrix_df,
    all_manifest_samples = all_manifest_samples,
    selected_samples = selected_samples,
    min_total_count = cfg$min_total_count,
    feature_type = feature_type
  )

  if (prepared$n_features < cfg$min_features_per_matrix) {
    skip_reason <- sprintf(
      "Insufficient features for %s: %d remaining after filtering (minimum required: %d). This matrix does not have enough signal for reliable DESeq2 analysis.",
      comparison, prepared$n_features, cfg$min_features_per_matrix
    )
    cat(sprintf("\n⚠️ ANALYSIS SKIPPED: %s - %s\n%s\n", comparison, feature_type, skip_reason))
    readr::write_lines(
      sprintf(
        "# ANALYSIS SKIPPED: %s - %s\n# Reason: Insufficient features after filtering\n# Features remaining: %d (minimum required: %d)",
        comparison, feature_type, prepared$n_features, cfg$min_features_per_matrix
      ),
      output_path
    )
    return(list(
      skipped = TRUE, skip_type = "insufficient_features", reason = skip_reason,
      feature_type = feature_type, display_name = category_display_name(feature_type),
      n_features = prepared$n_features
    ))
  }

  coldata <- selected_manifest %>%
    dplyr::select("sampleName", "groupName", "batch") %>%
    dplyr::mutate(group = factor(.data$groupName, levels = c(group2, group1))) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(coldata) <- coldata$sampleName

  dds <- DESeqDataSetFromMatrix(
    countData = prepared$counts_mat,
    colData = coldata,
    design = ~ group
  )
  dds <- estimateSizeFactors(dds, type = cfg$size_factor_type)

  fit_types_to_try <- c(cfg$fit_type, "local", "mean")
  fit_types_to_try <- unique(fit_types_to_try)
  dds_result <- NULL
  fit_errors <- list()

  for (fit_type in fit_types_to_try) {
    dds_result <- tryCatch(
      {
        DESeq(dds, fitType = fit_type, quiet = TRUE)
      },
      error = function(err) {
        fit_errors[[fit_type]] <<- conditionMessage(err)
        NULL
      }
    )

    if (!is.null(dds_result)) {
      if (fit_type != cfg$fit_type) {
        cat(sprintf(
          "\nNote: '%s' fit type failed. Successfully used '%s' as fallback.\n",
          cfg$fit_type, fit_type
        ))
      }
      dds <- dds_result
      break
    }
  }

  if (is.null(dds_result)) {
    n_group1 <- nrow(coldata[coldata$group == group1, ])
    n_group2 <- nrow(coldata[coldata$group == group2, ])
    failed_list <- paste(names(fit_errors), collapse = ", ")
    msg <- sprintf(
      paste(
        "\n=== DESeq2 DISPERSION ESTIMATION FAILURE ===",
        "Comparison: %s (%s vs %s)",
        "\nReplicates in this contrast:",
        "  %s: %d samples",
        "  %s: %d samples",
        "  (Recommended minimum: ≥3 per group)",
        "\nAll dispersion fitting methods failed: %s",
        "\nThis indicates a fundamental data quality issue:",
        "  • Sample groups are too similar or have too low variation",
        "  • Too few replicates per group (%d+%d total)",
        "  • Insufficient biological variation between samples",
        "  • Potential problems with the count matrix itself",
        "\nDebugging steps:",
        "  1. Verify this contrast is biologically meaningful",
        "  2. Check for data quality issues in the Tn5 matrices",
        "  3. Inspect actual peak counts for %s vs %s samples",
        "  4. Ensure samples are correctly labeled and distinct",
        "  5. Consider whether this comparison should be analyzed",
        "\nRecommendations:",
        "  • Add more replicates to this group",
        "  • Combine this contrast with others (if biologically justified)",
        "  • Skip this problematic contrast entirely",
        sep = "\n"
      ),
      comparison, group1, group2,
      group1, n_group1,
      group2, n_group2,
      failed_list,
      n_group1, n_group2,
      group1, group2
    )
    stop(msg, call. = FALSE)
  }

  contrast_vector <- c("group", group1, group2)
  res_raw <- results(
    dds,
    contrast = contrast_vector,
    alpha = cfg$alpha,
    pAdjustMethod = cfg$p_adjust_method,
    cooksCutoff = cfg$cooks_cutoff,
    independentFiltering = cfg$independent_filtering
  )

  shrinkage_error <- NULL
  res_shrunk <- tryCatch(
    {
      if (identical(cfg$shrink_type, "apeglm")) {
        lfcShrink(dds, coef = resultsNames(dds)[2], res = res_raw, type = "apeglm", quiet = TRUE)
      } else {
        lfcShrink(dds, contrast = contrast_vector, res = res_raw, type = cfg$shrink_type, quiet = TRUE)
      }
    },
    error = function(err) {
      shrinkage_error <<- conditionMessage(err)
      cat(sprintf(
        "\nNote: LFC shrinkage ('%s') failed for %s (%s): %s. Falling back to unshrunk estimates.\n",
        cfg$shrink_type, feature_type, comparison, conditionMessage(err)
      ))
      res_raw
    }
  )
  shrinkage_applied <- is.null(shrinkage_error)

  normalized_counts <- counts(dds, normalized = TRUE)
  group1_samples <- rownames(coldata)[coldata$group == group1]
  group2_samples <- rownames(coldata)[coldata$group == group2]

  result_tbl <- tibble::as_tibble(prepared$metadata_df) %>%
    dplyr::mutate(
      comparison = comparison,
      group1 = group1,
      group2 = group2,
      baseMean = res_raw$baseMean,
      log2FoldChange_raw = res_raw$log2FoldChange,
      log2FoldChange_shrunk = res_shrunk$log2FoldChange,
      shrinkage_applied = shrinkage_applied,
      lfcSE = res_raw$lfcSE,
      stat = res_raw$stat,
      pvalue = res_raw$pvalue,
      padj = res_raw$padj,
      # DESeq2's padj can be NA (e.g. independent filtering or Cook's outlier
      # removal) even when pvalue is present, silently dropping those features
      # from any padj-based significance call. adjusted_pvalue is a BH
      # correction recomputed directly from pvalue so significance can be
      # determined without relying on padj's NAs.
      adjusted_pvalue = p.adjust(res_raw$pvalue, method = "BH"),
      mean_group1 = rowMeans(normalized_counts[, group1_samples, drop = FALSE]),
      mean_group2 = rowMeans(normalized_counts[, group2_samples, drop = FALSE])
    ) %>%
    dplyr::select(
      "comparison",
      "group1",
      "group2",
      dplyr::everything()
    )

  readr::write_tsv(result_tbl, output_path)

  feature_label_col <- pick_feature_label_column(result_tbl)
  feature_labels <- safe_feature_labels(result_tbl, feature_label_col)

  significant_mask <- !is.na(result_tbl$adjusted_pvalue) &
    result_tbl$adjusted_pvalue < cfg$alpha &
    !is.na(result_tbl$log2FoldChange_shrunk) &
    abs(result_tbl$log2FoldChange_shrunk) >= cfg$lfc_threshold

  summary_list <- list(
    n_features_tested = nrow(result_tbl),
    n_significant = sum(significant_mask),
    n_up_group1 = sum(significant_mask & result_tbl$log2FoldChange_shrunk > 0, na.rm = TRUE),
    n_up_group2 = sum(significant_mask & result_tbl$log2FoldChange_shrunk < 0, na.rm = TRUE)
  )

  pca_data <- NULL
  percent_var <- NULL
  sample_distance <- NULL
  vsd <- tryCatch(vst(dds, blind = cfg$vst_blind), error = function(exc) NULL)
  if (!is.null(vsd)) {
    pca_raw <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
    percent_var <- round(100 * attr(pca_raw, "percentVar"), 2)
    pca_data <- as.data.frame(pca_raw)
    pca_data$sample <- rownames(pca_data)
    sample_distance <- as.matrix(dist(t(assay(vsd))))
  }

  list(
    skipped = FALSE,
    feature_type = feature_type,
    display_name = category_display_name(feature_type),
    matrix_path = matrix_path,
    output_path = output_path,
    feature_label_col = feature_label_col,
    feature_labels = feature_labels,
    results = result_tbl,
    summary = summary_list,
    shrinkage_applied = shrinkage_applied,
    shrinkage_error = shrinkage_error,
    pca_data = pca_data,
    percent_var = percent_var,
    sample_distance = sample_distance
  )
}

option_list <- list(
  make_option("--config-yaml", dest = "config_yaml", type = "character"),
  make_option("--manifest", dest = "manifest", type = "character"),
  make_option("--contrast-file", dest = "contrast_file", type = "character"),
  make_option("--comparison", dest = "comparison", type = "character"),
  make_option("--group1", dest = "group1", type = "character"),
  make_option("--group2", dest = "group2", type = "character"),
  make_option("--report-template", dest = "report_template", type = "character"),
  make_option("--report-output", dest = "report_output", type = "character"),
  make_option("--host-matrix", dest = "host_matrix", type = "character", default = ""),
  make_option("--host-output", dest = "host_output", type = "character", default = ""),
  make_option("--host-volcano", dest = "host_volcano", type = "character", default = ""),
  make_option("--virus-labels", dest = "virus_labels", type = "character", default = ""),
  make_option("--virus-matrices", dest = "virus_matrices", type = "character", default = ""),
  make_option("--virus-outputs", dest = "virus_outputs", type = "character", default = ""),
  make_option("--virus-volcanos", dest = "virus_volcanos", type = "character", default = "")
)

opts <- parse_args(OptionParser(option_list = option_list))

required_opts <- c(
  "config_yaml",
  "manifest",
  "contrast_file",
  "comparison",
  "group1",
  "group2",
  "report_template",
  "report_output"
)
for (opt_name in required_opts) {
  if (!nzchar(default_if_null(opts[[opt_name]], ""))) {
    stopf("Missing required option --%s", gsub("_", "-", opt_name))
  }
}

cfg <- read_config(opts$config_yaml)
validate_design_factors(default_if_null(cfg$deseq2$design_factors, c("group")))
deseq2_cfg <- build_deseq2_config(cfg)

manifest_data <- read_manifest(opts$manifest, opts$group1, opts$group2)
virus_labels <- split_delimited(opts$virus_labels)
virus_matrices <- split_delimited(opts$virus_matrices)
virus_outputs <- split_delimited(opts$virus_outputs)
virus_volcanos <- split_delimited(opts$virus_volcanos)

if (length(virus_labels) != length(virus_matrices) ||
    length(virus_labels) != length(virus_outputs) ||
    (length(virus_volcanos) > 0 && length(virus_volcanos) != length(virus_labels))) {
  stopf("Virus labels, matrix paths, output paths, and volcano paths must have the same length.")
}

available_category_keys <- c(
  if (nzchar(trimws(opts$host_matrix))) "host" else NULL,
  virus_labels
)
warn_unmatched_skip_features(deseq2_cfg$skip_features, available_category_keys)

dir.create(dirname(opts$report_output), recursive = TRUE, showWarnings = FALSE)

host_section <- NULL
if (nzchar(trimws(opts$host_matrix))) {
  if (!nzchar(trimws(opts$host_output))) {
    stopf("--host-output is required when --host-matrix is provided.")
  }
  dir.create(dirname(opts$host_output), recursive = TRUE, showWarnings = FALSE)
  host_section <- tryCatch(
    {
      run_deseq2_matrix(
        matrix_path = opts$host_matrix,
        output_path = opts$host_output,
        manifest = manifest_data$manifest,
        selected_manifest = manifest_data$selected,
        comparison = opts$comparison,
        group1 = opts$group1,
        group2 = opts$group2,
        feature_type = "host",
        cfg = deseq2_cfg
      )
    },
    error = function(err) {
      stopf("Host analysis failed for comparison '%s': %s", opts$comparison, conditionMessage(err))
    }
  )
  if (!is.null(host_section) && !isTRUE(host_section$skipped) && nzchar(trimws(opts$host_volcano))) {
    save_enhanced_volcano_png(
      result_tbl = host_section$results, feature_labels = host_section$feature_labels,
      png_path = opts$host_volcano, comparison = opts$comparison,
      group1 = opts$group1, group2 = opts$group2,
      alpha = deseq2_cfg$alpha, lfc_threshold = deseq2_cfg$lfc_threshold,
      label_top_n = deseq2_cfg$report_label_top_n
    )
  }
}

virus_sections <- list()
if (length(virus_labels) > 0) {
  for (idx in seq_along(virus_labels)) {
    dir.create(dirname(virus_outputs[[idx]]), recursive = TRUE, showWarnings = FALSE)
    virus_sections[[virus_labels[[idx]]]] <- tryCatch(
      {
        run_deseq2_matrix(
          matrix_path = virus_matrices[[idx]],
          output_path = virus_outputs[[idx]],
          manifest = manifest_data$manifest,
          selected_manifest = manifest_data$selected,
          comparison = opts$comparison,
          group1 = opts$group1,
          group2 = opts$group2,
          feature_type = virus_labels[[idx]],
          cfg = deseq2_cfg
        )
      },
      error = function(err) {
        stopf(
          "Analysis failed for %s in comparison '%s': %s",
          virus_labels[[idx]], opts$comparison, conditionMessage(err)
        )
      }
    )
    if (length(virus_volcanos) == length(virus_labels) &&
        !isTRUE(virus_sections[[virus_labels[[idx]]]]$skipped) &&
        nzchar(trimws(virus_volcanos[[idx]]))) {
      save_enhanced_volcano_png(
        result_tbl = virus_sections[[virus_labels[[idx]]]]$results,
        feature_labels = virus_sections[[virus_labels[[idx]]]]$feature_labels,
        png_path = virus_volcanos[[idx]], comparison = opts$comparison,
        group1 = opts$group1, group2 = opts$group2,
        alpha = deseq2_cfg$alpha, lfc_threshold = deseq2_cfg$lfc_threshold,
        label_top_n = deseq2_cfg$report_label_top_n
      )
    }
  }
}

context <- list(
  comparison = opts$comparison,
  group1 = opts$group1,
  group2 = opts$group2,
  contrast_file = opts$contrast_file,
  sample_table = manifest_data$selected %>%
    dplyr::select("sampleName", "groupName", "batch") %>%
    dplyr::rename(sample = "sampleName", group = "groupName"),
  config = deseq2_cfg,
  host = host_section,
  viruses = virus_sections
)

context_path <- tempfile(pattern = "deseq2_context_", fileext = ".rds")
on.exit(unlink(context_path), add = TRUE)
saveRDS(context, context_path)

rmarkdown::render(
  input = opts$report_template,
  output_file = basename(opts$report_output),
  output_dir = dirname(opts$report_output),
  params = list(context_path = context_path),
  quiet = TRUE,
  envir = new.env(parent = globalenv()),
  output_options = list(self_contained = deseq2_cfg$html_self_contained)
)
