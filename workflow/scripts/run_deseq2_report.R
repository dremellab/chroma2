#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
  library(rmarkdown)
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

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg) || !is.list(cfg)) stopf("Config file did not parse: %s", path)
  cfg
}

build_deseq2_config <- function(cfg) {
  dcfg <- default_if_null(cfg$deseq2, list())
  list(
    alpha = as.numeric(default_if_null(dcfg$alpha, 0.05)),
    lfc_threshold = as.numeric(default_if_null(dcfg$lfc_threshold, 1.0)),
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

option_list <- list(
  make_option("--config-yaml", dest = "config_yaml", type = "character"),
  make_option("--manifest", dest = "manifest", type = "character"),
  make_option("--contrast-file", dest = "contrast_file", type = "character"),
  make_option("--tsv-dir", dest = "tsv_dir", type = "character"),
  make_option("--report-template", dest = "report_template", type = "character"),
  make_option("--output", dest = "output", type = "character")
)

opts <- parse_args(OptionParser(option_list = option_list))

required_opts <- c("config_yaml", "manifest", "contrast_file", "tsv_dir", "report_template", "output")
for (opt_name in required_opts) {
  if (!nzchar(default_if_null(opts[[opt_name]], ""))) stopf("Missing required option --%s", gsub("_", "-", opt_name))
}

cfg <- read_config(opts$config_yaml)
deseq2_cfg <- build_deseq2_config(cfg)

contrasts_df <- suppressMessages(readr::read_tsv(opts$contrast_file, show_col_types = FALSE, progress = FALSE))
if (list(contrasts_df$columns) != list(c("group1", "group2"))) {
  stopf("contrasts_file must have exactly two headers: group1 and group2")
}

contrasts_list <- list()
for (idx in seq_len(nrow(contrasts_df))) {
  group1 <- contrasts_df[[idx, "group1"]]
  group2 <- contrasts_df[[idx, "group2"]]
  comparison <- paste(group1, "vs", group2, sep = "_")
  contrasts_list[[comparison]] <- list(group1 = group1, group2 = group2)
}

extract_skip_info <- function(tsv_file) {
  lines <- readLines(tsv_file, n = 10)
  skip_info <- list(skipped = FALSE, reason = NA_character_, skip_type = NA_character_)

  for (line in lines) {
    if (grepl("^# ANALYSIS SKIPPED", line)) {
      skip_info$skipped <- TRUE
    }
    if (grepl("^# Reason:", line)) {
      skip_info$reason <- trimws(sub("^# Reason:", "", line))
    }
    if (grepl("^# Skip list:", line)) {
      skip_info$skip_type <- "skip_list"
    }
    if (grepl("^# Reason: Insufficient features", line)) {
      skip_info$skip_type <- "insufficient_features"
    }
  }

  skip_info
}

tsv_files <- list.files(opts$tsv_dir, pattern = ".*\\.deseq2\\.tsv$", full.names = TRUE, recursive = TRUE)
if (length(tsv_files) == 0) stopf("No DESeq2 TSV files found in %s", opts$tsv_dir)

result_sections <- list()
for (tsv_file in tsv_files) {
  result_tbl <- suppressMessages(readr::read_tsv(tsv_file, show_col_types = FALSE, progress = FALSE, comment = "#"))

  comparison <- NA_character_
  group1 <- NA_character_
  group2 <- NA_character_

  if (nrow(result_tbl) > 0) {
    comparison <- result_tbl$comparison[[1]]
    group1 <- result_tbl$group1[[1]]
    group2 <- result_tbl$group2[[1]]
    is_skipped <- FALSE
  } else {
    is_skipped <- TRUE
    lines <- readLines(tsv_file, n = 1)
    if (length(lines) > 0 && grepl("^# ANALYSIS SKIPPED:", lines[1])) {
      parts <- strsplit(lines[1], " - ")[[1]]
      if (length(parts) >= 2) {
        comparison_part <- trimws(sub("^# ANALYSIS SKIPPED: ", "", parts[1]))
        comparison <- comparison_part
      }
    }
    if (is.na(comparison)) {
      next
    }
    group1 <- NA_character_
    group2 <- NA_character_
  }

  if (is.na(comparison)) next

  if (!comparison %in% names(result_sections)) {
    result_sections[[comparison]] <- list(
      group1 = group1, group2 = group2, matrices = list()
    )
  } else if (!is.na(group1) && !is.na(group2)) {
    result_sections[[comparison]]$group1 <- group1
    result_sections[[comparison]]$group2 <- group2
  }

  matrix_label <- basename(tools::file_path_sans_ext(tsv_file))
  matrix_label <- gsub(paste0("^", comparison, "\\."), "", matrix_label)
  matrix_label <- gsub("\\.deseq2$", "", matrix_label)

  if (is_skipped) {
    skip_info <- extract_skip_info(tsv_file)
    result_sections[[comparison]]$matrices[[matrix_label]] <- list(
      skipped = TRUE,
      skip_info = skip_info,
      tsv_path = tsv_file
    )
  } else {
    result_sections[[comparison]]$matrices[[matrix_label]] <- list(
      skipped = FALSE,
      results = result_tbl,
      tsv_path = tsv_file
    )
  }
}

context <- list(
  contrasts = contrasts_list,
  result_sections = result_sections,
  config = deseq2_cfg,
  matrix_labels = unique(unlist(lapply(result_sections, function(x) names(x$matrices))))
)

context_path <- tempfile(pattern = "deseq2_context_", fileext = ".rds")
on.exit(unlink(context_path), add = TRUE)
saveRDS(context, context_path)

rmarkdown::render(
  input = opts$report_template,
  output_file = basename(opts$output),
  output_dir = dirname(opts$output),
  params = list(context_path = context_path),
  quiet = TRUE,
  envir = new.env(parent = globalenv()),
  output_options = list(self_contained = deseq2_cfg$html_self_contained)
)

cat(sprintf("Report saved to %s\n", opts$output))
