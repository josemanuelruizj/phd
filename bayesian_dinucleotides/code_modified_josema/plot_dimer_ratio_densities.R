#!/usr/bin/env Rscript

parse_args <- function(args) {
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(
      "Usage:\n",
      "  Rscript plot_dimer_ratio_densities.R --input_dir results_all_genomes_1kb --methylation_table dnmt.tsv --output_dir dimer_density_plots\n\n",
      "Options:\n",
      "  --input_dir           Directory containing *_dinucleotides_by_window.tsv files. Required unless --density_table is used.\n",
      "  --methylation_table   TSV with columns species, taxonomic_group_from_figure, methylation_CpG_OE. Optional.\n",
      "  --output_dir          Output directory. Default: dimer_density_plots\n",
      "  --density_table       Existing density TSV to plot without reading raw dimer tables. Optional.\n",
      "  --ratio_col           Ratio column to use. Default: ratio_observed_expected\n",
      "  --x_min               Minimum x value for density grid. Default: 0\n",
      "  --x_max               Maximum x value for density grid. Default: 2.2\n",
      "  --grid_n              Number of density grid points. Default: 512\n",
      "  --bw                  Density bandwidth passed to density(). Default: nrd0\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  opts <- list(
    input_dir = NULL,
    methylation_table = NULL,
    output_dir = "dimer_density_plots",
    density_table = NULL,
    ratio_col = "ratio_observed_expected",
    x_min = 0,
    x_max = 2.2,
    grid_n = 512L,
    bw = "nrd0"
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

  if (is.null(opts$density_table) && is.null(opts$input_dir)) {
    stop("Use --input_dir to compute densities or --density_table to reuse existing densities", call. = FALSE)
  }
  if (!is.null(opts$input_dir) && !dir.exists(opts$input_dir)) {
    stop("Input directory does not exist: ", opts$input_dir, call. = FALSE)
  }
  if (!is.null(opts$density_table) && !file.exists(opts$density_table)) {
    stop("Density table does not exist: ", opts$density_table, call. = FALSE)
  }
  if (!is.null(opts$methylation_table) && !file.exists(opts$methylation_table)) {
    stop("Methylation table does not exist: ", opts$methylation_table, call. = FALSE)
  }

  opts$x_min <- as.numeric(opts$x_min)
  opts$x_max <- as.numeric(opts$x_max)
  opts$grid_n <- as.integer(opts$grid_n)
  if (opts$x_min >= opts$x_max) stop("--x_min must be lower than --x_max", call. = FALSE)
  if (opts$grid_n < 64L) stop("--grid_n must be at least 64", call. = FALSE)

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
  files <- list.files(
    input_dir,
    pattern = "_dinucleotides_by_window\\.tsv(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  sort(files)
}

compute_one_density <- function(values, x_min, x_max, grid_n, bw) {
  values <- values[is.finite(values)]
  if (length(values) < 2L || length(unique(values)) < 2L) {
    return(data.frame(x = numeric(), density = numeric()))
  }
  den <- stats::density(values, from = x_min, to = x_max, n = grid_n, bw = bw, na.rm = TRUE)
  data.frame(x = den$x, density = den$y)
}

compute_density_table <- function(opts) {
  files <- list_dimer_files(opts$input_dir)
  if (length(files) == 0L) {
    stop("No *_dinucleotides_by_window.tsv files found in: ", opts$input_dir, call. = FALSE)
  }

  parts <- vector("list", length(files))
  for (i in seq_along(files)) {
    path <- files[[i]]
    message("Reading ", basename(path), " [", i, "/", length(files), "]")
    tab <- read_tsv(path, select = c("species", "dimer", opts$ratio_col))
    missing <- setdiff(c("dimer", opts$ratio_col), names(tab))
    if (length(missing) > 0L) {
      stop("File ", path, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    species_id <- if ("species" %in% names(tab) && any(!is.na(tab$species))) {
      as.character(tab$species[which(!is.na(tab$species))[1L]])
    } else {
      species_from_file(path)
    }
    tab$dimer_group <- dimer_group(tab$dimer)
    tab[[opts$ratio_col]] <- as.numeric(tab[[opts$ratio_col]])

    group_parts <- lapply(split(tab[[opts$ratio_col]], tab$dimer_group), function(x) {
      den <- compute_one_density(x, opts$x_min, opts$x_max, opts$grid_n, opts$bw)
      if (nrow(den) == 0L) return(NULL)
      den$n_values <- sum(is.finite(x))
      den
    })
    group_parts <- group_parts[!vapply(group_parts, is.null, logical(1L))]
    if (length(group_parts) == 0L) next

    out <- do.call(rbind, Map(function(group_name, den) {
      den$species_id <- species_id
      den$species_display <- species_display_name(species_id)
      den$dimer_group <- group_name
      den$source_file <- basename(path)
      den
    }, names(group_parts), group_parts))
    rownames(out) <- NULL
    parts[[i]] <- out
    rm(tab)
    gc(verbose = FALSE)
  }

  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) stop("No density values could be computed", call. = FALSE)
  out <- do.call(rbind, parts)
  out$dimer_group <- factor(out$dimer_group, levels = dimer_group_levels)
  out <- out[!is.na(out$dimer_group), , drop = FALSE]
  out[order(out$dimer_group, out$species_display, out$x), , drop = FALSE]
}

add_metadata <- function(density, methylation_table) {
  density$species_key <- species_key(density$species_display)
  density$taxonomic_group_from_figure <- "Unknown"
  density$methylation_CpG_OE <- "unknown"

  if (is.null(methylation_table)) return(density)

  meta <- read_tsv(methylation_table)
  if (!"species" %in% names(meta)) stop("Methylation table must contain a species column", call. = FALSE)
  meta$species_key <- species_key(meta$species)

  keep <- intersect(c("species_key", "taxonomic_group_from_figure", "methylation_CpG_OE"), names(meta))
  meta <- unique(meta[, keep, drop = FALSE])
  out <- merge(density, meta, by = "species_key", all.x = TRUE, suffixes = c("", ".meta"))

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

make_plot <- function(density, color_col, title, output_prefix, output_dir, x_min, x_max) {
  p <- ggplot2::ggplot(
    density,
    ggplot2::aes(x = x, y = density, group = species_id, color = .data[[color_col]])
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    ggplot2::geom_line(alpha = 0.85, linewidth = 0.35) +
    ggplot2::facet_grid(dimer_group ~ ., scales = "free_y") +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max), expand = FALSE) +
    ggplot2::labs(
      title = title,
      x = "Observed/Expected Dinucleotide Ratio",
      y = "Density",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey80", color = "grey30"),
      strip.text.y = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  pdf_path <- file.path(output_dir, paste0(output_prefix, ".pdf"))
  png_path <- file.path(output_dir, paste0(output_prefix, ".png"))
  ggplot2::ggsave(pdf_path, p, width = 11, height = 12, limitsize = FALSE)
  ggplot2::ggsave(png_path, p, width = 11, height = 12, dpi = 300, limitsize = FALSE)
  message("Wrote: ", pdf_path)
  message("Wrote: ", png_path)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  require_pkg("ggplot2")
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

  density_path <- if (!is.null(opts$density_table)) {
    opts$density_table
  } else {
    file.path(opts$output_dir, "dimer_ratio_density_values.tsv")
  }

  if (!is.null(opts$density_table)) {
    density <- read_tsv(density_path)
  } else {
    density <- compute_density_table(opts)
    write_tsv(density, density_path)
  }

  density <- add_metadata(density, opts$methylation_table)
  density$dimer_group <- factor(as.character(density$dimer_group), levels = dimer_group_levels)
  density <- density[!is.na(density$dimer_group), , drop = FALSE]
  density$species_display <- factor(density$species_display, levels = sort(unique(density$species_display)))
  density$taxonomic_group_from_figure <- factor(density$taxonomic_group_from_figure)
  density$methylation_CpG_OE <- factor(density$methylation_CpG_OE)

  annotated_path <- file.path(opts$output_dir, "dimer_ratio_density_values_with_metadata.tsv")
  write_tsv(density, annotated_path)

  make_plot(
    density,
    "species_display",
    "Distribution of Dimer Ratio Values Across Species",
    "dimer_ratio_density_by_species",
    opts$output_dir,
    opts$x_min,
    opts$x_max
  )
  make_plot(
    density,
    "taxonomic_group_from_figure",
    "Distribution of Dimer Ratio Values Across Clades",
    "dimer_ratio_density_by_clade",
    opts$output_dir,
    opts$x_min,
    opts$x_max
  )
  make_plot(
    density,
    "methylation_CpG_OE",
    "Distribution of Dimer Ratio Values by CpG Methylation Status",
    "dimer_ratio_density_by_methylation",
    opts$output_dir,
    opts$x_min,
    opts$x_max
  )
}

main()
