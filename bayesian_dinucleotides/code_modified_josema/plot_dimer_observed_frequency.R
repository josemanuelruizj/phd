#!/usr/bin/env Rscript

parse_args <- function(args) {
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(
      "Usage:\n",
      "  Rscript plot_dimer_observed_frequency.R --input_dir results_all_genomes_1kb --methylation_table dnmt.tsv --output_dir dimer_observed_frequency_plots\n\n",
      "Options:\n",
      "  --input_dir           Directory containing *_dinucleotides_by_window.tsv files. Required unless --observed_table is used.\n",
      "  --methylation_table   TSV with species, taxonomic_group_from_figure, methylation_CpG_OE. Optional.\n",
      "  --output_dir          Output directory. Default: dimer_observed_frequency_plots\n",
      "  --observed_table      Existing observed-frequency TSV to plot without reading raw tables. Optional.\n",
      "  --frequency_col       Observed-frequency column to use. Default: observed_frequency\n",
      "  --x_max               Optional max x value for density plots. Default: auto per whole plot\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  opts <- list(
    input_dir = NULL,
    methylation_table = NULL,
    output_dir = "dimer_observed_frequency_plots",
    observed_table = NULL,
    frequency_col = "observed_frequency",
    x_max = NA_real_
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key, call. = FALSE)
    name <- sub("^--", "", key)
    if (!name %in% names(opts)) stop("Unknown option: --", name, call. = FALSE)
    if (i == length(args)) stop("Missing value for option: --", name, call. = FALSE)
    opts[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }

  if (is.null(opts$observed_table) && is.null(opts$input_dir)) {
    stop("Use --input_dir to build the observed table or --observed_table to reuse an existing table", call. = FALSE)
  }
  if (!is.null(opts$input_dir) && !dir.exists(opts$input_dir)) {
    stop("Input directory does not exist: ", opts$input_dir, call. = FALSE)
  }
  if (!is.null(opts$observed_table) && !file.exists(opts$observed_table)) {
    stop("Observed table does not exist: ", opts$observed_table, call. = FALSE)
  }
  if (!is.null(opts$methylation_table) && !file.exists(opts$methylation_table)) {
    stop("Methylation table does not exist: ", opts$methylation_table, call. = FALSE)
  }
  opts$x_max <- as.numeric(opts$x_max)
  opts
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

read_tsv <- function(path, select = NULL) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, sep = "\t", select = select, showProgress = FALSE)))
  }
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  tab <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(select)) tab <- tab[, intersect(select, names(tab)), drop = FALSE]
  tab
}

write_tsv <- function(tab, path) {
  write.table(tab, path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Wrote: ", path)
}

species_from_file <- function(path) {
  x <- basename(path)
  sub("_dinucleotides_by_window\\.tsv(\\.gz)?$", "", x, ignore.case = TRUE)
}

species_display_name <- function(species_id) {
  x <- sub("__.*$", "", species_id)
  gsub("_", " ", x, fixed = TRUE)
}

species_key <- function(x) {
  tolower(gsub("[_[:space:]]+", " ", trimws(x)))
}

revcomp_dimer <- function(x) {
  chars <- strsplit(x, "", fixed = TRUE)[[1L]]
  paste(rev(chartr("ACGT", "TGCA", chars)), collapse = "")
}

dimer_group <- function(dimer) {
  dimer <- toupper(dimer)
  rc <- vapply(dimer, revcomp_dimer, character(1L))
  pair <- ifelse(dimer <= rc, paste(dimer, rc, sep = "/"), paste(rc, dimer, sep = "/"))
  pair[dimer == rc] <- dimer[dimer == rc]
  pair
}

dimer_group_levels <- c("CG", "GC", "CC/GG", "AG/CT", "GA/TC", "CA/TG", "AC/GT", "AA/TT", "TA", "AT")

list_dimer_files <- function(input_dir) {
  sort(list.files(
    input_dir,
    pattern = "_dinucleotides_by_window\\.tsv(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  ))
}

build_observed_table <- function(opts) {
  files <- list_dimer_files(opts$input_dir)
  if (length(files) == 0L) {
    stop("No *_dinucleotides_by_window.tsv files found in: ", opts$input_dir, call. = FALSE)
  }

  parts <- vector("list", length(files))
  for (i in seq_along(files)) {
    path <- files[[i]]
    message("Reading ", basename(path), " [", i, "/", length(files), "]")
    tab <- read_tsv(path, select = c("species", "seq_id", "window_index", "start", "end", "dimer", opts$frequency_col, "observed_count"))
    missing <- setdiff(c("dimer", opts$frequency_col), names(tab))
    if (length(missing) > 0L) {
      stop("File ", path, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    species_id <- if ("species" %in% names(tab) && any(!is.na(tab$species))) {
      as.character(tab$species[which(!is.na(tab$species))[1L]])
    } else {
      species_from_file(path)
    }

    tab$dimer_group <- dimer_group(tab$dimer)
    tab$dimer_observed_frequency <- as.numeric(tab[[opts$frequency_col]])
    tab$species_id <- species_id
    tab$species_display <- species_display_name(species_id)
    tab$source_file <- basename(path)

    keep <- intersect(
      c("species_id", "species_display", "seq_id", "window_index", "start", "end", "dimer", "dimer_group", "observed_count", "dimer_observed_frequency", "source_file"),
      names(tab)
    )
    parts[[i]] <- tab[, keep, drop = FALSE]
    rm(tab)
    gc(verbose = FALSE)
  }

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) stop("No dimer rows found", call. = FALSE)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out$dimer_group <- factor(out$dimer_group, levels = dimer_group_levels)
  out <- out[!is.na(out$dimer_group), , drop = FALSE]
  out
}

add_metadata <- function(tab, methylation_table) {
  tab$species_key <- species_key(tab$species_display)
  tab$taxonomic_group_from_figure <- "Unknown"
  tab$methylation_CpG_OE <- "unknown"

  if (is.null(methylation_table)) return(tab)

  meta <- read_tsv(methylation_table)
  if (!"species" %in% names(meta)) stop("Methylation table must contain a species column", call. = FALSE)
  meta$species_key <- species_key(meta$species)
  keep <- intersect(c("species_key", "taxonomic_group_from_figure", "methylation_CpG_OE"), names(meta))
  meta <- unique(meta[, keep, drop = FALSE])
  out <- merge(tab, meta, by = "species_key", all.x = TRUE, suffixes = c("", ".meta"))

  if ("taxonomic_group_from_figure.meta" %in% names(out)) {
    out$taxonomic_group_from_figure <- ifelse(
      is.na(out$taxonomic_group_from_figure.meta) | out$taxonomic_group_from_figure.meta == "",
      "Unknown",
      out$taxonomic_group_from_figure.meta
    )
    out$taxonomic_group_from_figure.meta <- NULL
  }
  if ("methylation_CpG_OE.meta" %in% names(out)) {
    out$methylation_CpG_OE <- ifelse(
      is.na(out$methylation_CpG_OE.meta) | out$methylation_CpG_OE.meta == "",
      "unknown",
      out$methylation_CpG_OE.meta
    )
    out$methylation_CpG_OE.meta <- NULL
  }
  out
}

summarise_species_dimer <- function(tab) {
  key <- interaction(tab$species_id, tab$dimer_group, drop = TRUE, lex.order = TRUE)
  parts <- lapply(split(tab, key), function(x) {
    data.frame(
      species_id = x$species_id[[1L]],
      species_display = x$species_display[[1L]],
      dimer_group = as.character(x$dimer_group[[1L]]),
      taxonomic_group_from_figure = x$taxonomic_group_from_figure[[1L]],
      methylation_CpG_OE = x$methylation_CpG_OE[[1L]],
      n_values = sum(is.finite(x$dimer_observed_frequency)),
      mean_observed_frequency = mean(x$dimer_observed_frequency, na.rm = TRUE),
      median_observed_frequency = stats::median(x$dimer_observed_frequency, na.rm = TRUE),
      q25_observed_frequency = stats::quantile(x$dimer_observed_frequency, 0.25, na.rm = TRUE, names = FALSE),
      q75_observed_frequency = stats::quantile(x$dimer_observed_frequency, 0.75, na.rm = TRUE, names = FALSE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out$dimer_group <- factor(out$dimer_group, levels = dimer_group_levels)
  out[order(out$dimer_group, out$methylation_CpG_OE, out$species_display), , drop = FALSE]
}

save_plot <- function(plot, output_dir, prefix, width = 11, height = 12) {
  pdf_path <- file.path(output_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(output_dir, paste0(prefix, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  message("Wrote: ", pdf_path)
  message("Wrote: ", png_path)
}

plot_density <- function(tab, color_col, title, x_max) {
  ggplot2::ggplot(tab, ggplot2::aes(x = dimer_observed_frequency, color = .data[[color_col]], group = species_id)) +
    ggplot2::geom_density(linewidth = 0.35, alpha = 0.85, na.rm = TRUE) +
    ggplot2::facet_grid(dimer_group ~ ., scales = "free_y") +
    ggplot2::coord_cartesian(xlim = c(0, x_max), expand = FALSE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = title,
      x = "Observed dinucleotide frequency per window",
      y = "Density",
      color = NULL
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey80", color = "grey30"),
      strip.text.y = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_violin <- function(tab) {
  ggplot2::ggplot(tab, ggplot2::aes(x = methylation_CpG_OE, y = dimer_observed_frequency, fill = methylation_CpG_OE)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.65, na.rm = TRUE) +
    ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.2, alpha = 0.85, na.rm = TRUE) +
    ggplot2::facet_wrap(~ dimer_group, scales = "free_y", ncol = 2) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = "Observed Dinucleotide Frequency by CpG Methylation Status",
      x = "CpG methylation status",
      y = "Observed dinucleotide frequency per window",
      fill = NULL
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), legend.position = "none")
}

plot_species_summary <- function(summary_tab) {
  ggplot2::ggplot(
    summary_tab,
    ggplot2::aes(x = species_display, y = median_observed_frequency, color = methylation_CpG_OE)
  ) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ dimer_group, scales = "free_x", ncol = 2) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(
      title = "Median Observed Dinucleotide Frequency per Species",
      x = NULL,
      y = "Median observed frequency per window",
      color = "CpG methylation"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  require_pkg("ggplot2")
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

  observed_path <- if (!is.null(opts$observed_table)) opts$observed_table else file.path(opts$output_dir, "dimer_observed_frequency_by_window.tsv")
  if (!is.null(opts$observed_table)) {
    observed <- read_tsv(observed_path)
  } else {
    observed <- build_observed_table(opts)
    write_tsv(observed, observed_path)
  }

  observed <- add_metadata(observed, opts$methylation_table)
  observed$dimer_observed_frequency <- as.numeric(observed$dimer_observed_frequency)
  observed$dimer_group <- factor(as.character(observed$dimer_group), levels = dimer_group_levels)
  observed$methylation_CpG_OE <- factor(observed$methylation_CpG_OE)
  observed$taxonomic_group_from_figure <- factor(observed$taxonomic_group_from_figure)
  observed$species_display <- factor(observed$species_display, levels = sort(unique(observed$species_display)))
  observed <- observed[!is.na(observed$dimer_group), , drop = FALSE]

  annotated_path <- file.path(opts$output_dir, "dimer_observed_frequency_by_window_with_metadata.tsv")
  write_tsv(observed, annotated_path)

  summary_tab <- summarise_species_dimer(observed)
  write_tsv(summary_tab, file.path(opts$output_dir, "dimer_observed_frequency_by_species.tsv"))

  x_max <- opts$x_max
  if (!is.finite(x_max)) {
    x_max <- stats::quantile(observed$dimer_observed_frequency, 0.995, na.rm = TRUE, names = FALSE)
    x_max <- max(x_max, 0.001, na.rm = TRUE)
  }

  save_plot(
    plot_density(observed, "species_display", "Observed Dinucleotide Frequency Across Species", x_max),
    opts$output_dir,
    "dimer_observed_frequency_density_by_species",
    width = 11,
    height = 12
  )
  save_plot(
    plot_density(observed, "taxonomic_group_from_figure", "Observed Dinucleotide Frequency Across Clades", x_max),
    opts$output_dir,
    "dimer_observed_frequency_density_by_clade",
    width = 10,
    height = 12
  )
  save_plot(
    plot_density(observed, "methylation_CpG_OE", "Observed Dinucleotide Frequency by CpG Methylation Status", x_max),
    opts$output_dir,
    "dimer_observed_frequency_density_by_methylation",
    width = 9,
    height = 12
  )
  save_plot(
    plot_violin(observed),
    opts$output_dir,
    "dimer_observed_frequency_violin_by_methylation",
    width = 10,
    height = 12
  )
  save_plot(
    plot_species_summary(summary_tab),
    opts$output_dir,
    "dimer_observed_frequency_species_medians",
    width = 12,
    height = 13
  )
}

main()
