#!/usr/bin/env Rscript

parse_args <- function(args) {
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(
      "Usage:\n",
      "  Rscript annotate_functional_summaries_chipseeker.R --summary_dir summaries --annotation_dir annotations --output_dir functional_annotation\n\n",
      "Options:\n",
      "  --summary_dir       Directory with functional summary TSV files. Required.\n",
      "  --annotation_dir    Directory with GTF/GFF/GFF3/TxDb files. Required unless --annotation_file is used.\n",
      "  --annotation_file   One annotation file to use for one species. Optional.\n",
      "  --output_dir        Output directory. Default: functional_annotation\n",
      "  --fasta_dir         Optional FASTA directory used only to infer species names.\n",
      "  --species           Optional species name. Defaults to all species detected in --summary_dir.\n",
      "  --tss_upstream      Bases upstream of TSS for promoter annotation. Default: 3000\n",
      "  --tss_downstream    Bases downstream of TSS for promoter annotation. Default: 3000\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  opts <- list(
    summary_dir = NULL,
    annotation_dir = NULL,
    annotation_file = NULL,
    output_dir = "functional_annotation",
    fasta_dir = NULL,
    species = NULL,
    tss_upstream = 3000L,
    tss_downstream = 3000L
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

  if (is.null(opts$summary_dir)) stop("Required option missing: --summary_dir", call. = FALSE)
  if (!dir.exists(opts$summary_dir)) stop("Summary directory does not exist: ", opts$summary_dir, call. = FALSE)
  if (is.null(opts$annotation_file) && is.null(opts$annotation_dir)) {
    stop("Use --annotation_dir or --annotation_file", call. = FALSE)
  }
  if (!is.null(opts$annotation_dir) && !dir.exists(opts$annotation_dir)) {
    stop("Annotation directory does not exist: ", opts$annotation_dir, call. = FALSE)
  }
  if (!is.null(opts$annotation_file) && !file.exists(opts$annotation_file)) {
    stop("Annotation file does not exist: ", opts$annotation_file, call. = FALSE)
  }
  if (!is.null(opts$fasta_dir) && !dir.exists(opts$fasta_dir)) {
    stop("FASTA directory does not exist: ", opts$fasta_dir, call. = FALSE)
  }

  opts$tss_upstream <- as.integer(opts$tss_upstream)
  opts$tss_downstream <- as.integer(opts$tss_downstream)
  if (opts$tss_upstream < 0L || opts$tss_downstream < 0L) {
    stop("TSS distances must be >= 0", call. = FALSE)
  }

  opts
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

strip_fasta_ext <- function(path) {
  x <- basename(path)
  x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.(fasta|fa|fna|fas)$", "", x, ignore.case = TRUE)
  x
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

summary_label <- function(path, species) {
  x <- basename(path)
  x <- sub("\\.tsv$", "", x)
  sub(paste0("^", escape_regex(species), "_"), "", x)
}

species_from_summary_files <- function(summary_dir) {
  files <- list.files(
    summary_dir,
    pattern = "_functional_.*\\.tsv$",
    full.names = FALSE
  )
  files <- files[!grepl("_annotation_|_combined_annotation|_summary_by_|_genes_by_", files)]
  unique(sub("_functional_.*$", "", files))
}

species_from_fasta_dir <- function(fasta_dir) {
  if (is.null(fasta_dir)) return(character())
  files <- list.files(
    fasta_dir,
    pattern = "\\.(fa|fasta|fna|fas)(\\.gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  unique(vapply(files, strip_fasta_ext, character(1L)))
}

list_summary_files <- function(summary_dir, species) {
  files <- list.files(
    summary_dir,
    pattern = paste0("^", escape_regex(species), "_functional_.*\\.tsv$"),
    full.names = TRUE
  )
  files[!grepl("_annotation_|_combined_annotation|_summary_by_|_genes_by_", basename(files))]
}

annotation_candidates <- function(annotation_dir) {
  list.files(
    annotation_dir,
    pattern = "\\.(gtf|gff|gff3|sqlite|sqlite3|rds)(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

find_annotation_file <- function(species, annotation_dir) {
  candidates <- annotation_candidates(annotation_dir)
  if (length(candidates) == 0L) {
    stop("No GTF/GFF/GFF3/TxDb annotation files found in: ", annotation_dir, call. = FALSE)
  }

  base <- basename(candidates)
  species_re <- escape_regex(species)
  exactish <- grepl(paste0("^", species_re, "([._-]|$)"), base, ignore.case = TRUE)
  contains <- grepl(species_re, base, ignore.case = TRUE)

  matches <- candidates[exactish]
  if (length(matches) == 0L) matches <- candidates[contains]
  if (length(matches) == 0L) {
    stop(
      "No annotation file matched species '", species, "'. ",
      "Rename the annotation with the FASTA/species prefix or pass --annotation_file.",
      call. = FALSE
    )
  }
  if (length(matches) > 1L) {
    message("Multiple annotation files matched ", species, "; using: ", matches[[1L]])
  }
  matches[[1L]]
}

make_txdb <- function(annotation_file) {
  lower <- tolower(annotation_file)
  if (grepl("\\.rds$", lower)) {
    txdb <- readRDS(annotation_file)
    if (!inherits(txdb, "TxDb")) stop("--annotation_file RDS is not a TxDb object", call. = FALSE)
    return(txdb)
  }
  if (grepl("\\.(sqlite|sqlite3)$", lower)) {
    return(GenomicFeatures::loadDb(annotation_file))
  }
  if (grepl("\\.(gtf|gff|gff3)(\\.gz)?$", lower)) {
    return(GenomicFeatures::makeTxDbFromGFF(annotation_file))
  }
  stop("Unsupported annotation format: ", annotation_file, call. = FALSE)
}

make_ranges_from_summary <- function(tab, source_file) {
  required <- c("seq_id", "start", "end", "functional_group")
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0L) {
    stop("Summary table ", source_file, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  tab$region_id <- paste(tab$seq_id, tab$start, tab$end, tab$functional_group, seq_len(nrow(tab)), sep = "_")
  ranges <- GenomicRanges::GRanges(
    seqnames = tab$seq_id,
    ranges = IRanges::IRanges(start = as.integer(tab$start), end = as.integer(tab$end)),
    strand = "*"
  )
  S4Vectors::mcols(ranges) <- S4Vectors::DataFrame(tab)
  names(ranges) <- tab$region_id
  ranges
}

write_tsv <- function(table, path) {
  write.table(table, path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Wrote: ", path)
}

safe_annotation_group <- function(annotation) {
  out <- as.character(annotation)
  out[is.na(out) | out == ""] <- "Unannotated"
  sub(" \\(.*$", "", out)
}

count_percent <- function(tab, group_cols) {
  if (nrow(tab) == 0L) return(data.frame())
  key <- interaction(tab[group_cols], drop = TRUE, lex.order = TRUE)
  counts <- as.data.frame(table(key), stringsAsFactors = FALSE)
  names(counts) <- c("key", "count")
  parts <- do.call(rbind, strsplit(counts$key, ".", fixed = TRUE))
  out <- as.data.frame(parts, stringsAsFactors = FALSE)
  names(out) <- group_cols
  out$count <- counts$count

  total_key <- interaction(out[setdiff(group_cols, "annotation_group")], drop = TRUE, lex.order = TRUE)
  totals <- ave(out$count, total_key, FUN = sum)
  out$percentage <- out$count / totals * 100
  out
}

plot_outputs <- function(combined, output_dir, species) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(combined) == 0L) return(invisible(NULL))

  by_group <- count_percent(combined, c("functional_group", "annotation_group"))
  if (nrow(by_group) > 0L) {
    p1 <- ggplot2::ggplot(by_group, ggplot2::aes(x = functional_group, y = percentage, fill = annotation_group)) +
      ggplot2::geom_col() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Genomic annotation by CAN-CAN group:", species),
        x = "Functional CAN-CAN group",
        y = "Percentage"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file.path(output_dir, paste0(species, "_annotation_by_functional_group.pdf")), p1, width = 12, height = 8)
  }

  if ("distanceToTSS" %in% names(combined)) {
    p2 <- ggplot2::ggplot(combined, ggplot2::aes(x = functional_group, y = distanceToTSS, fill = functional_group)) +
      ggplot2::geom_boxplot(outlier.size = 0.4) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Distance to TSS by CAN-CAN group:", species),
        x = "Functional CAN-CAN group",
        y = "Distance to TSS (bp)"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
    ggplot2::ggsave(file.path(output_dir, paste0(species, "_distance_to_TSS_by_functional_group.pdf")), p2, width = 12, height = 8)
  }
}

process_species <- function(species, opts) {
  message("========================================")
  message("Processing species: ", species)

  summary_files <- list_summary_files(opts$summary_dir, species)
  if (length(summary_files) == 0L) {
    warning("No summary tables found for species: ", species, immediate. = TRUE)
    return(invisible(NULL))
  }

  annotation_file <- if (!is.null(opts$annotation_file)) {
    opts$annotation_file
  } else {
    find_annotation_file(species, opts$annotation_dir)
  }
  message("Annotation file: ", annotation_file)
  txdb <- make_txdb(annotation_file)

  species_dir <- file.path(opts$output_dir, species)
  dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)

  combined <- data.frame()
  for (summary_file in summary_files) {
    label <- summary_label(summary_file, species)
    message("Annotating summary table: ", basename(summary_file))
    tab <- read.delim(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(tab) == 0L) {
      warning("Skipping empty summary table: ", summary_file, immediate. = TRUE)
      next
    }

    peaks <- make_ranges_from_summary(tab, summary_file)
    anno <- ChIPseeker::annotatePeak(
      peaks,
      TxDb = txdb,
      tssRegion = c(-opts$tss_upstream, opts$tss_downstream),
      verbose = FALSE
    )
    anno_df <- as.data.frame(anno)
    anno_df$species <- species
    anno_df$summary_table <- basename(summary_file)
    anno_df$summary_label <- label
    anno_df$annotation_group <- safe_annotation_group(anno_df$annotation)

    out_file <- file.path(species_dir, paste0(species, "_", label, "_chipseeker_annotation.tsv"))
    write_tsv(anno_df, out_file)
    combined <- rbind(combined, anno_df)
  }

  if (nrow(combined) == 0L) return(invisible(NULL))

  combined_file <- file.path(species_dir, paste0(species, "_all_summary_tables_chipseeker_annotation.tsv"))
  write_tsv(combined, combined_file)

  by_functional_group <- count_percent(combined, c("summary_label", "functional_group", "annotation_group"))
  write_tsv(
    by_functional_group,
    file.path(species_dir, paste0(species, "_summary_by_functional_group.tsv"))
  )

  by_summary_table <- count_percent(combined, c("summary_label", "annotation_group"))
  write_tsv(
    by_summary_table,
    file.path(species_dir, paste0(species, "_summary_by_summary_table.tsv"))
  )

  gene_col <- intersect(c("SYMBOL", "gene_name", "geneId", "gene_id"), names(combined))
  if (length(gene_col) > 0L) {
    gene_col <- gene_col[[1L]]
    genes <- combined[!is.na(combined[[gene_col]]) & combined[[gene_col]] != "", , drop = FALSE]
    keep <- intersect(
      c("species", "summary_label", "functional_group", gene_col, "annotation", "distanceToTSS"),
      names(genes)
    )
    genes <- unique(genes[, keep, drop = FALSE])
    write_tsv(genes, file.path(species_dir, paste0(species, "_genes_by_functional_group.tsv")))
  }

  plot_outputs(combined, species_dir, species)
  invisible(combined_file)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  require_pkg("ChIPseeker")
  require_pkg("GenomicFeatures")
  require_pkg("GenomicRanges")
  require_pkg("IRanges")
  require_pkg("S4Vectors")

  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

  species <- if (!is.null(opts$species)) {
    opts$species
  } else {
    from_fasta <- species_from_fasta_dir(opts$fasta_dir)
    from_summaries <- species_from_summary_files(opts$summary_dir)
    if (length(from_fasta) > 0L) intersect(from_fasta, from_summaries) else from_summaries
  }

  if (length(species) == 0L) {
    stop("No species detected. Check --summary_dir or pass --species.", call. = FALSE)
  }

  message("Species to process: ", paste(species, collapse = ", "))
  for (sp in species) process_species(sp, opts)
  message("Functional annotation finished. Output: ", opts$output_dir)
}

main()
