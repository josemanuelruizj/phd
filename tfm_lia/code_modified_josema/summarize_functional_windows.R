#!/usr/bin/env Rscript

parse_args <- function(args) {
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(
      "Usage:\n",
      "  Rscript summarize_functional_windows.R --functional_table species_bhlh_CAN_CAN_functional_by_window.tsv --output_dir summaries\n\n",
      "Options:\n",
      "  --functional_table       Functional CAN-CAN table produced by genome_window_analysis.R. Required.\n",
      "  --output_dir             Output directory for summary tables. Default: summaries\n",
      "  --species                Output prefix. Defaults to functional table basename.\n",
      "  --top_n                  Top N windows per functional group and direction. Default: 100\n",
      "  --top_percent            Top percent windows per functional group and direction. Default: 1\n",
      "  --log2_threshold         Absolute pseudocount log2 threshold. Default: 1\n",
      "  --pvalue_threshold       Adjusted p-value threshold. Default: 0.05\n",
      "  --min_expected_zero      Minimum expected count for zero-observed table. Default: 1\n",
      "  --pseudocount            Pseudocount for stable log2 ratios. Default: 0.5\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  opts <- list(
    functional_table = NULL,
    output_dir = "summaries",
    species = NULL,
    top_n = 100L,
    top_percent = 1,
    log2_threshold = 1,
    pvalue_threshold = 0.05,
    min_expected_zero = 1,
    pseudocount = 0.5
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (!name %in% names(opts)) {
      stop("Unknown option: --", name, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for option: --", name, call. = FALSE)
    }
    opts[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }

  if (is.null(opts$functional_table)) {
    stop("Required option missing: --functional_table", call. = FALSE)
  }
  if (!file.exists(opts$functional_table)) {
    stop("Functional table does not exist: ", opts$functional_table, call. = FALSE)
  }
  if (is.null(opts$species)) {
    base <- basename(opts$functional_table)
    opts$species <- sub("_bhlh_CAN_CAN_functional_by_window\\.tsv$", "", base)
    opts$species <- sub("\\.tsv$", "", opts$species)
  }

  opts$top_n <- as.integer(opts$top_n)
  opts$top_percent <- as.numeric(opts$top_percent)
  opts$log2_threshold <- as.numeric(opts$log2_threshold)
  opts$pvalue_threshold <- as.numeric(opts$pvalue_threshold)
  opts$min_expected_zero <- as.numeric(opts$min_expected_zero)
  opts$pseudocount <- as.numeric(opts$pseudocount)

  if (opts$top_n < 1L) stop("--top_n must be at least 1", call. = FALSE)
  if (opts$top_percent <= 0 || opts$top_percent > 100) {
    stop("--top_percent must be > 0 and <= 100", call. = FALSE)
  }
  if (opts$log2_threshold < 0) stop("--log2_threshold must be 0 or higher", call. = FALSE)
  if (opts$pvalue_threshold <= 0 || opts$pvalue_threshold > 1) {
    stop("--pvalue_threshold must be > 0 and <= 1", call. = FALSE)
  }
  if (opts$min_expected_zero < 0) stop("--min_expected_zero must be 0 or higher", call. = FALSE)
  if (opts$pseudocount <= 0) stop("--pseudocount must be > 0", call. = FALSE)

  opts
}

expected_column <- function(table) {
  if ("expected_count_markov" %in% names(table)) {
    "expected_count_markov"
  } else if ("expected_count" %in% names(table)) {
    "expected_count"
  } else {
    stop("Input table must contain expected_count_markov or expected_count", call. = FALSE)
  }
}

ensure_metrics <- function(table, pseudocount) {
  expected_col <- expected_column(table)
  expected <- table[[expected_col]]
  observed <- table$observed_count

  table$ratio_observed_expected <- ifelse(expected > 0, observed / expected, NA_real_)
  table$log2_ratio_observed_expected <- ifelse(
    is.na(table$ratio_observed_expected),
    NA_real_,
    log2(table$ratio_observed_expected)
  )
  table$log2_ratio_pseudocount <- log2((observed + pseudocount) / (expected + pseudocount))
  table$depletion_score <- expected - observed
  table$enrichment_score <- observed - expected
  table$is_zero_observed <- observed == 0 & expected > 0
  table$p_under <- ifelse(expected > 0, ppois(observed, lambda = expected), NA_real_)
  table$p_over <- ifelse(expected > 0, ppois(observed - 1, lambda = expected, lower.tail = FALSE), NA_real_)
  table$padj_under <- p.adjust(table$p_under, method = "BH")
  table$padj_over <- p.adjust(table$p_over, method = "BH")

  table$padj_under_by_group <- ave(
    table$p_under,
    table$functional_group,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  table$padj_over_by_group <- ave(
    table$p_over,
    table$functional_group,
    FUN = function(x) p.adjust(x, method = "BH")
  )

  table
}

write_tsv <- function(table, path) {
  write.table(table, path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Wrote: ", path)
}

add_direction <- function(table, direction) {
  if (nrow(table) == 0L) {
    return(table)
  }
  table$representation <- direction
  table
}

top_by_group <- function(table, n_or_percent, percent = FALSE) {
  parts <- lapply(split(table, table$functional_group), function(group_table) {
    n_keep <- if (percent) {
      max(1L, ceiling(nrow(group_table) * n_or_percent / 100))
    } else {
      n_or_percent
    }

    over <- group_table[group_table$log2_ratio_pseudocount > 0, , drop = FALSE]
    under <- group_table[group_table$log2_ratio_pseudocount < 0, , drop = FALSE]
    over <- over[order(-over$log2_ratio_pseudocount, over$padj_over), , drop = FALSE]
    under <- under[order(under$log2_ratio_pseudocount, under$padj_under), , drop = FALSE]

    rbind(
      add_direction(head(over, n_keep), "over"),
      add_direction(head(under, n_keep), "under")
    )
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

depletion_top_by_group <- function(table, top_n) {
  parts <- lapply(split(table, table$functional_group), function(group_table) {
    under <- group_table[group_table$depletion_score > 0, , drop = FALSE]
    under <- under[order(-under$depletion_score, under$padj_under), , drop = FALSE]
    add_direction(head(under, top_n), "under")
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

  functional <- read.delim(opts$functional_table, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("functional_group", "observed_count")
  missing <- setdiff(required, names(functional))
  if (length(missing) > 0L) {
    stop("Functional table is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  functional <- ensure_metrics(functional, opts$pseudocount)
  functional <- functional[!is.na(functional$log2_ratio_pseudocount), , drop = FALSE]
  percent_label <- gsub(".", "p", as.character(opts$top_percent), fixed = TRUE)
  threshold_label <- gsub(".", "p", as.character(opts$log2_threshold), fixed = TRUE)
  pvalue_label <- gsub(".", "p", as.character(opts$pvalue_threshold), fixed = TRUE)

  write_tsv(
    top_by_group(functional, opts$top_n, percent = FALSE),
    file.path(opts$output_dir, paste0(opts$species, "_functional_topN_by_log2_pseudocount.tsv"))
  )

  write_tsv(
    depletion_top_by_group(functional, opts$top_n),
    file.path(opts$output_dir, paste0(opts$species, "_functional_topN_by_depletion_score.tsv"))
  )

  write_tsv(
    top_by_group(functional, opts$top_percent, percent = TRUE),
    file.path(opts$output_dir, paste0(opts$species, "_functional_top", percent_label, "percent_by_log2_pseudocount.tsv"))
  )

  log2_threshold <- functional[abs(functional$log2_ratio_pseudocount) >= opts$log2_threshold, , drop = FALSE]
  log2_threshold$representation <- ifelse(log2_threshold$log2_ratio_pseudocount >= 0, "over", "under")
  log2_threshold <- log2_threshold[order(
    log2_threshold$functional_group,
    -abs(log2_threshold$log2_ratio_pseudocount)
  ), , drop = FALSE]
  write_tsv(
    log2_threshold,
    file.path(opts$output_dir, paste0(opts$species, "_functional_log2_pseudocount_threshold_", threshold_label, ".tsv"))
  )

  pvalue_significant <- rbind(
    add_direction(functional[functional$padj_over <= opts$pvalue_threshold, , drop = FALSE], "over"),
    add_direction(functional[functional$padj_under <= opts$pvalue_threshold, , drop = FALSE], "under")
  )
  pvalue_significant <- pvalue_significant[order(
    pvalue_significant$representation,
    pvalue_significant$functional_group,
    pmin(pvalue_significant$padj_over, pvalue_significant$padj_under, na.rm = TRUE)
  ), , drop = FALSE]
  write_tsv(
    pvalue_significant,
    file.path(opts$output_dir, paste0(opts$species, "_functional_padj_significant_", pvalue_label, ".tsv"))
  )

  expected_col <- expected_column(functional)
  zero_observed <- functional[
    functional$is_zero_observed & functional[[expected_col]] >= opts$min_expected_zero,
    ,
    drop = FALSE
  ]
  zero_observed <- zero_observed[order(-zero_observed[[expected_col]], zero_observed$padj_under), , drop = FALSE]
  zero_observed$representation <- "under"
  write_tsv(
    zero_observed,
    file.path(opts$output_dir, paste0(opts$species, "_functional_zero_observed_high_expected.tsv"))
  )
}

main()
