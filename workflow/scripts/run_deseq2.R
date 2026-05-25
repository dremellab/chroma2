#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(EnhancedVolcano)
})

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

default_if_null <- function(value, default) {
  if (is.null(value)) return(default)
  value
}

as_bool <- function(value, name) {
  if (is.logical(value) && length(value) == 1 && !is.na(value)) return(value)
  if (is.numeric(value) && length(value) == 1 && !is.na(value)) return(value != 0)
  if (is.character(value) && length(value) == 1) {
    lowered <- tolower(trimws(value))
    if (lowered %in% c("true", "1", "yes")) return(TRUE)
    if (lowered %in% c("false", "0", "no")) return(FALSE)
  }
  stopf("%s must be a boolean value. Found: %s", name, paste(value, collapse = ","))
}

pick_feature_label_column <- function(df) {
  candidates <- c("gene_name", "gene_id", "bin_id", "chrom")
  for (candidate in candidates) {
    if (candidate %in% colnames(df)) {
      values <- as.character(df[[candidate]])
      if (any(nzchar(values))) return(candidate)
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
  padj_safe <- ifelse(is.na(result_tbl$padj), 1, result_tbl$padj)

  sig_mask <- !is.na(result_tbl$padj) &
    result_tbl$padj < alpha &
    !is.na(result_tbl$log2FoldChange_shrunk) &
    abs(result_tbl$log2FoldChange_shrunk) >= lfc_threshold

  select_lab <- rep("", nrow(result_tbl))
  if (sum(sig_mask) > 0) {
    sig_idx <- which(sig_mask)
    ord <- order(result_tbl$padj[sig_idx], -abs(result_tbl$log2FoldChange_shrunk[sig_idx]))
    top_idx <- sig_idx[ord][seq_len(min(label_top_n, length(sig_idx)))]
    select_lab[top_idx] <- feature_labels[top_idx]
  }

  p <- EnhancedVolcano::EnhancedVolcano(
    toptable = data.frame(log2FC = result_tbl$log2FoldChange_shrunk,
                          padj = padj_safe, stringsAsFactors = FALSE),
    lab = select_lab, x = "log2FC", y = "padj",
    title = sprintf("%s  |  %s vs %s", comparison, group1, group2),
    subtitle = bquote(alpha == .(alpha) ~ "| LFC ==" ~ .(lfc_threshold)),
    pCutoff = alpha, FCcutoff = lfc_threshold,
    xlab = bquote(Shrunken ~ log[2] ~ fold ~ change),
    ylab = bquote(-log[10] ~ p[adj]),
    legendPosition = "right", drawConnectors = TRUE,
    widthConnectors = 0.5, colConnectors = "grey30",
    pointSize = 1.5, labSize = 3.0, max.overlaps = Inf
  )

  ggplot2::ggsave(filename = png_path, plot = p,
                  width = 10, height = 8, dpi = 150, units = "in")
  invisible(png_path)
}

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg) || !is.list(cfg)) stopf("Config file did not parse: %s", path)
  cfg
}

build_deseq2_config <- function(cfg) {
  dcfg <- default_if_null(cfg$deseq2, list())
  skip_features_raw <- default_if_null(dcfg$skip_features, list())
  if (is.character(skip_features_raw)) {
    skip_features <- tolower(as.character(skip_features_raw))
  } else if (is.list(skip_features_raw)) {
    skip_features <- tolower(unlist(skip_features_raw, use.names = FALSE))
  } else {
    skip_features <- character(0)
  }
  list(
    alpha = as.numeric(default_if_null(dcfg$alpha, 0.05)),
    lfc_threshold = as.numeric(default_if_null(dcfg$lfc_threshold, 1.0)),
    min_total_count = as.integer(default_if_null(dcfg$min_total_count, 1)),
    min_features_per_matrix = as.integer(default_if_null(dcfg$min_features_per_matrix, 50)),
    skip_features = skip_features,
    size_factor_type = as.character(default_if_null(dcfg$size_factor_type, "poscounts")),
    fit_type = as.character(default_if_null(dcfg$fit_type, "parametric")),
    shrink_type = as.character(default_if_null(dcfg$shrink_type, "apeglm")),
    p_adjust_method = as.character(default_if_null(dcfg$p_adjust_method, "BH")),
    cooks_cutoff = as_bool(default_if_null(dcfg$cooks_cutoff, TRUE), "deseq2.cooks_cutoff"),
    independent_filtering = as_bool(default_if_null(dcfg$independent_filtering, TRUE), "deseq2.independent_filtering"),
    vst_blind = as_bool(default_if_null(dcfg$vst_blind, TRUE), "deseq2.vst_blind"),
    report_top_n = as.integer(default_if_null(dcfg$report_top_n, 50)),
    report_label_top_n = as.integer(default_if_null(dcfg$report_label_top_n, 15)),
    report_max_table_rows = as.integer(default_if_null(dcfg$report_max_table_rows, 50000)),
    html_self_contained = as_bool(default_if_null(dcfg$html_self_contained, TRUE), "deseq2.html_self_contained")
  )
}

read_manifest <- function(path, group1, group2) {
  manifest <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  required_cols <- c("sampleName", "groupName", "batch")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) stopf("Manifest missing columns: %s", paste(missing_cols, collapse = ", "))

  selected <- manifest %>%
    dplyr::filter(.data$groupName %in% c(group1, group2)) %>%
    dplyr::mutate(
      sampleName = as.character(.data$sampleName),
      groupName = as.character(.data$groupName),
      batch = as.character(.data$batch)
    )

  if (nrow(selected) == 0) stopf("No manifest samples matched groups '%s' and '%s'.", group1, group2)
  if (!all(c(group1, group2) %in% selected$groupName)) {
    stopf("Manifest selection did not retain both groups '%s' and '%s'.", group1, group2)
  }

  list(manifest = manifest, selected = selected)
}

prepare_counts <- function(matrix_df, all_manifest_samples, selected_samples, min_total_count, feature_type) {
  metadata_cols <- setdiff(colnames(matrix_df), all_manifest_samples)
  missing_samples <- setdiff(selected_samples, colnames(matrix_df))
  if (length(missing_samples) > 0) stopf("Matrix missing columns: %s", paste(missing_samples, collapse = ", "))

  counts_df <- matrix_df[, selected_samples, drop = FALSE]
  counts_mat <- as.matrix(counts_df)
  storage.mode(counts_mat) <- "double"
  if (any(is.na(counts_mat))) stopf("Matrix contains NA counts.")
  if (any(counts_mat < 0)) stopf("Matrix contains negative counts.")
  if (any(abs(counts_mat - round(counts_mat)) > 1e-8)) stopf("Matrix contains non-integer counts.")
  counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"

  total_features_before <- nrow(counts_mat)
  row_totals <- rowSums(counts_mat)
  keep <- row_totals >= min_total_count
  features_filtered <- sum(!keep)
  features_remaining <- sum(keep)

  if (!any(keep)) stopf("No features passed min_total_count=%d.", min_total_count)

  if (features_filtered > 0) {
    message("[INFO] Filtering features with total counts < ", min_total_count)
    message("[INFO]   Total features before filtering: ", total_features_before)
    message("[INFO]   Features filtered: ", features_filtered)
    message("[INFO]   Features remaining: ", features_remaining)

    filtered_features <- which(!keep)
    if (length(filtered_features) <= 20) {
      filter_info <- paste0(
        "Feature: ", rownames(matrix_df)[filtered_features],
        " (total counts: ", row_totals[filtered_features], ")"
      )
      for (msg in filter_info) {
        message("[DEBUG]     ", msg)
      }
    } else {
      count_dist <- table(cut(row_totals[!keep], breaks = c(-Inf, 1, 5, 10, Inf)))
      message("[DEBUG]     Distribution of filtered feature counts:")
      for (i in seq_along(count_dist)) {
        message("[DEBUG]       ", names(count_dist)[i], ": ", count_dist[i], " features")
      }
    }
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

  list(metadata_df = metadata_df, counts_mat = counts_mat, n_features = features_remaining)
}

run_deseq2_matrix <- function(
  matrix_path, output_path, volcano_path, manifest, selected_manifest,
  comparison, group1, group2, feature_type, cfg
) {
  feature_type_lower <- tolower(feature_type)

  if (feature_type_lower %in% cfg$skip_features) {
    skip_reason <- sprintf(
      "Feature type '%s' is in the skip list (configured: %s)",
      feature_type, paste(cfg$skip_features, collapse = ", ")
    )
    message("")
    message("================================================================================")
    message("⏭️  ANALYSIS SKIPPED: ", comparison, " - ", feature_type)
    message("================================================================================")
    message(skip_reason)
    message("")

    skip_header <- sprintf(
      "# ANALYSIS SKIPPED: %s - %s\n# Reason: Feature type in skip list\n# Skip list: %s",
      comparison, feature_type, paste(cfg$skip_features, collapse = ", ")
    )
    readr::write_lines(skip_header, output_path)
    readr::write_lines(skip_header, volcano_path)

    return(list(skipped = TRUE, reason = skip_reason, skip_type = "skip_list", n_features = NA))
  }

  matrix_df <- suppressMessages(readr::read_tsv(matrix_path, show_col_types = FALSE, progress = FALSE))
  if (nrow(matrix_df) == 0) stopf("Matrix has no rows: %s", matrix_path)

  all_manifest_samples <- as.character(manifest$sampleName)
  selected_samples <- as.character(selected_manifest$sampleName)
  prepared <- prepare_counts(
    matrix_df = matrix_df, all_manifest_samples = all_manifest_samples,
    selected_samples = selected_samples, min_total_count = cfg$min_total_count,
    feature_type = feature_type
  )

  if (prepared$n_features < cfg$min_features_per_matrix) {
    warning_msg <- sprintf(
      "⚠️ ANALYSIS SKIPPED: Insufficient features for %s\n    Features remaining after filtering: %d\n    Minimum required: %d\n    This matrix does not have enough features to reliably estimate dispersion in DESeq2.",
      comparison, prepared$n_features, cfg$min_features_per_matrix
    )
    message("")
    message("================================================================================")
    message(warning_msg)
    message("================================================================================")
    message("")

    skip_header <- sprintf(
      "# ANALYSIS SKIPPED: %s - %s\n# Reason: Insufficient features after filtering\n# Features remaining: %d (minimum required: %d)\n# This matrix does not have enough signal for reliable DESeq2 analysis.",
      comparison, feature_type, prepared$n_features, cfg$min_features_per_matrix
    )
    readr::write_lines(skip_header, output_path)
    readr::write_lines(skip_header, volcano_path)

    return(list(skipped = TRUE, reason = warning_msg, skip_type = "insufficient_features", n_features = prepared$n_features))
  }

  coldata <- selected_manifest %>%
    dplyr::select("sampleName", "groupName", "batch") %>%
    dplyr::mutate(group = factor(.data$groupName, levels = c(group2, group1))) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(coldata) <- coldata$sampleName

  dds <- DESeqDataSetFromMatrix(
    countData = prepared$counts_mat, colData = coldata, design = ~ group
  )
  dds <- estimateSizeFactors(dds, type = cfg$size_factor_type)

  dds <- tryCatch(
    {
      DESeq(dds, fitType = cfg$fit_type, quiet = TRUE)
    },
    error = function(e) {
      message("Standard dispersion estimation failed, using gene-wise estimates: ", e$message)
      dds <<- estimateDispersionsGeneEst(dds)
      dispersions(dds) <<- mcols(dds)$dispGeneEst
      dds <<- nbinomWaldTest(dds)
      dds
    }
  )

  contrast_vector <- c("group", group1, group2)
  res_raw <- results(
    dds, contrast = contrast_vector, alpha = cfg$alpha,
    pAdjustMethod = cfg$p_adjust_method, cooksCutoff = cfg$cooks_cutoff,
    independentFiltering = cfg$independent_filtering
  )

  res_shrunk <- res_raw
  normalized_counts <- counts(dds, normalized = TRUE)
  group1_samples <- rownames(coldata)[coldata$group == group1]
  group2_samples <- rownames(coldata)[coldata$group == group2]

  result_tbl <- tibble::as_tibble(prepared$metadata_df) %>%
    dplyr::mutate(
      comparison = comparison, group1 = group1, group2 = group2,
      baseMean = res_raw$baseMean,
      log2FoldChange_raw = res_raw$log2FoldChange,
      log2FoldChange_shrunk = res_shrunk$log2FoldChange,
      lfcSE = res_raw$lfcSE, stat = res_raw$stat,
      pvalue = res_raw$pvalue, padj = res_raw$padj,
      mean_group1 = rowMeans(normalized_counts[, group1_samples, drop = FALSE]),
      mean_group2 = rowMeans(normalized_counts[, group2_samples, drop = FALSE])
    ) %>%
    dplyr::select("comparison", "group1", "group2", dplyr::everything())

  readr::write_tsv(result_tbl, output_path)

  feature_label_col <- pick_feature_label_column(result_tbl)
  feature_labels <- safe_feature_labels(result_tbl, feature_label_col)

  save_enhanced_volcano_png(
    result_tbl = result_tbl, feature_labels = feature_labels,
    png_path = volcano_path, comparison = comparison,
    group1 = group1, group2 = group2,
    alpha = cfg$alpha, lfc_threshold = cfg$lfc_threshold,
    label_top_n = cfg$report_label_top_n
  )

  vsd <- tryCatch(vst(dds, blind = cfg$vst_blind), error = function(exc) NULL)
  pca_data <- NULL
  percent_var <- NULL
  if (!is.null(vsd)) {
    pca_raw <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
    percent_var <- round(100 * attr(pca_raw, "percentVar"), 2)
    pca_data <- as.data.frame(pca_raw)
    pca_data$sample <- rownames(pca_data)
  }

  list(
    skipped = FALSE,
    feature_type = feature_type, matrix_path = matrix_path,
    output_path = output_path, volcano_path = volcano_path,
    feature_label_col = feature_label_col, feature_labels = feature_labels,
    results = result_tbl, pca_data = pca_data, percent_var = percent_var
  )
}

option_list <- list(
  make_option("--config-yaml", dest = "config_yaml", type = "character"),
  make_option("--manifest", dest = "manifest", type = "character"),
  make_option("--matrix", dest = "matrix", type = "character"),
  make_option("--output", dest = "output", type = "character"),
  make_option("--volcano", dest = "volcano", type = "character"),
  make_option("--comparison", dest = "comparison", type = "character"),
  make_option("--group1", dest = "group1", type = "character"),
  make_option("--group2", dest = "group2", type = "character"),
  make_option("--matrix-label", dest = "matrix_label", type = "character"),
  make_option("--feature-type", dest = "feature_type", type = "character")
)

opts <- parse_args(OptionParser(option_list = option_list))

required_opts <- c("config_yaml", "manifest", "matrix", "output", "volcano", "comparison", "group1", "group2", "matrix_label", "feature_type")
for (opt_name in required_opts) {
  if (!nzchar(default_if_null(opts[[opt_name]], ""))) stopf("Missing required option --%s", gsub("_", "-", opt_name))
}

cfg <- read_config(opts$config_yaml)
deseq2_cfg <- build_deseq2_config(cfg)
manifest_data <- read_manifest(opts$manifest, opts$group1, opts$group2)

message("[INFO] Running DESeq2 analysis")
message("[INFO]   Comparison: ", opts$comparison)
message("[INFO]   Matrix label: ", opts$matrix_label)
message("[INFO]   Groups: ", opts$group1, " vs ", opts$group2)
message("[INFO]   Min total count per feature: ", deseq2_cfg$min_total_count)
message("[INFO]   Min features required: ", deseq2_cfg$min_features_per_matrix)
if (length(deseq2_cfg$skip_features) > 0) {
  message("[INFO]   Skip features: ", paste(deseq2_cfg$skip_features, collapse = ", "))
}

result <- run_deseq2_matrix(
  matrix_path = opts$matrix, output_path = opts$output, volcano_path = opts$volcano,
  manifest = manifest_data$manifest, selected_manifest = manifest_data$selected,
  comparison = opts$comparison, group1 = opts$group1, group2 = opts$group2,
  feature_type = opts$feature_type, cfg = deseq2_cfg
)

if (result$skipped) {
  skip_type_msg <- if (result$skip_type == "skip_list") {
    "[SKIP] User-configured skip list"
  } else if (result$skip_type == "insufficient_features") {
    "[SKIP] Insufficient features after filtering"
  } else {
    "[SKIP] Other reason"
  }
  message(skip_type_msg, ": ", opts$comparison, " - ", opts$matrix_label)
  quit(save = "no", status = 0)
} else {
  message("[INFO] DESeq2 analysis completed successfully")
}
