#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
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

split_delimited <- function(value, delimiter = ";;;") {
  value <- trimws(default_if_null(value, ""))
  if (identical(value, "")) {
    return(character(0))
  }
  pieces <- strsplit(value, delimiter, fixed = TRUE)[[1]]
  trimws(pieces[nzchar(trimws(pieces))])
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

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg) || !is.list(cfg)) {
    stopf("Config file did not parse into a mapping: %s", path)
  }
  cfg
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
    counts_mat = counts_mat
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
  dds <- DESeq(dds, fitType = cfg$fit_type, quiet = TRUE)

  contrast_vector <- c("group", group1, group2)
  res_raw <- results(
    dds,
    contrast = contrast_vector,
    alpha = cfg$alpha,
    pAdjustMethod = cfg$p_adjust_method,
    cooksCutoff = cfg$cooks_cutoff,
    independentFiltering = cfg$independent_filtering
  )

  res_shrunk <- res_raw

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
      lfcSE = res_raw$lfcSE,
      stat = res_raw$stat,
      pvalue = res_raw$pvalue,
      padj = res_raw$padj,
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

  significant_mask <- !is.na(result_tbl$padj) &
    result_tbl$padj < cfg$alpha &
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
    feature_type = feature_type,
    matrix_path = matrix_path,
    output_path = output_path,
    feature_label_col = feature_label_col,
    feature_labels = feature_labels,
    results = result_tbl,
    summary = summary_list,
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
  make_option("--virus-labels", dest = "virus_labels", type = "character", default = ""),
  make_option("--virus-matrices", dest = "virus_matrices", type = "character", default = ""),
  make_option("--virus-outputs", dest = "virus_outputs", type = "character", default = "")
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

if (length(virus_labels) != length(virus_matrices) ||
    length(virus_labels) != length(virus_outputs)) {
  stopf("Virus labels, matrix paths, and output paths must have the same length.")
}

dir.create(dirname(opts$report_output), recursive = TRUE, showWarnings = FALSE)

host_section <- NULL
if (nzchar(trimws(opts$host_matrix))) {
  if (!nzchar(trimws(opts$host_output))) {
    stopf("--host-output is required when --host-matrix is provided.")
  }
  dir.create(dirname(opts$host_output), recursive = TRUE, showWarnings = FALSE)
  host_section <- run_deseq2_matrix(
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
}

virus_sections <- list()
if (length(virus_labels) > 0) {
  for (idx in seq_along(virus_labels)) {
    dir.create(dirname(virus_outputs[[idx]]), recursive = TRUE, showWarnings = FALSE)
    virus_sections[[virus_labels[[idx]]]] <- run_deseq2_matrix(
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
