#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(parallel)
})

parse_args <- function(args) {
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(
      "Usage:\n",
      "  Rscript genome_window_analysis.R --fasta genome.fa.gz --species human --output_dir results --window_size 50000 --cores 8\n\n",
      "Options:\n",
      "  --fasta              Complete reference genome FASTA (.fa, .fasta, optionally .gz). Required.\n",
      "  --species            Species/sample label. Defaults to FASTA basename.\n",
      "  --output_dir         Output directory. Default: results\n",
      "  --window_size        Window size in bp. Default: 50000\n",
      "  --step_size          Window step in bp. Default: window_size\n",
      "  --cores              CPUs to use inside one node. Default: 1\n",
      "  --include_revcomp    TRUE/FALSE, analyse reverse complement too. Default: TRUE\n",
      "  --min_acgt_fraction  Drop windows below this A/C/G/T fraction. Default: 0\n",
      "  --top_n              Top N over/underrepresented windows per functional group. Default: 100\n",
      "  --log2_threshold     Absolute log2 ratio threshold for functional regions. Default: 1\n",
      "  --top_percent        Top percent over/underrepresented windows per functional group. Default: 1\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }

  opts <- list(
    fasta = NULL,
    species = NULL,
    output_dir = "results",
    window_size = 50000L,
    step_size = NA_integer_,
    cores = 1L,
    include_revcomp = TRUE,
    min_acgt_fraction = 0,
    top_n = 100L,
    log2_threshold = 1,
    top_percent = 1
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
    value <- args[[i + 1L]]
    opts[[name]] <- value
    i <- i + 2L
  }

  if (is.null(opts$fasta)) {
    stop("Required option missing: --fasta", call. = FALSE)
  }

  if (is.null(opts$species)) {
    base <- basename(opts$fasta)
    opts$species <- sub("\\.(fa|fasta)(\\.gz)?$", "", base, ignore.case = TRUE)
  }

  opts$window_size <- as.integer(opts$window_size)
  opts$cores <- as.integer(opts$cores)
  opts$include_revcomp <- tolower(opts$include_revcomp) %in% c("true", "t", "1", "yes", "y")
  opts$min_acgt_fraction <- as.numeric(opts$min_acgt_fraction)
  opts$top_n <- as.integer(opts$top_n)
  opts$log2_threshold <- as.numeric(opts$log2_threshold)
  opts$top_percent <- as.numeric(opts$top_percent)

  if (is.na(as.integer(opts$step_size))) {
    opts$step_size <- opts$window_size
  } else {
    opts$step_size <- as.integer(opts$step_size)
  }

  if (opts$window_size < 6L) {
    stop("--window_size must be at least 6", call. = FALSE)
  }
  if (opts$step_size < 1L) {
    stop("--step_size must be at least 1", call. = FALSE)
  }
  if (opts$cores < 1L) {
    stop("--cores must be at least 1", call. = FALSE)
  }
  if (opts$min_acgt_fraction < 0 || opts$min_acgt_fraction > 1) {
    stop("--min_acgt_fraction must be between 0 and 1", call. = FALSE)
  }
  if (opts$top_n < 1L) {
    stop("--top_n must be at least 1", call. = FALSE)
  }
  if (opts$log2_threshold < 0) {
    stop("--log2_threshold must be 0 or higher", call. = FALSE)
  }
  if (opts$top_percent <= 0 || opts$top_percent > 100) {
    stop("--top_percent must be > 0 and <= 100", call. = FALSE)
  }

  opts
}

log2_or_na <- function(x) {
  ifelse(is.na(x), NA_real_, log2(x))
}

open_fasta <- function(path) {
  if (!file.exists(path)) {
    stop("FASTA file does not exist: ", path, call. = FALSE)
  }
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

sanitize_seq_name <- function(header) {
  sub("[[:space:]].*$", "", sub("^>", "", header))
}

reverse_complement <- function(seq) {
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  paste(rev(chartr("ACGTNacgtn", "TGCANtgcan", chars)), collapse = "")
}

chunk_windows <- function(seq_len, window_size, step_size) {
  starts <- seq.int(1L, seq_len, by = step_size)
  ends <- pmin(starts + window_size - 1L, seq_len)
  data.frame(window_index = seq_along(starts), start = starts, end = ends)
}

count_kmers <- function(seq, k, alphabet_regex) {
  seq_len <- nchar(seq)
  if (seq_len < k) {
    return(integer(0))
  }
  kmers <- substring(seq, 1:(seq_len - k + 1L), k:seq_len)
  kmers <- kmers[grepl(alphabet_regex, kmers)]
  table(kmers)
}

dinucleotide_table_for_window <- function(window_seq, meta, include_revcomp) {
  seq <- toupper(window_seq)
  analysis_seq <- if (include_revcomp) {
    paste0(seq, "NNNNNN", reverse_complement(seq))
  } else {
    seq
  }

  dimers <- as.vector(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0))
  bases <- strsplit(analysis_seq, "", fixed = TRUE)[[1]]
  acgt_bases <- bases[bases %in% c("A", "C", "G", "T")]
  base_counts <- table(factor(acgt_bases, levels = c("A", "C", "G", "T")))
  base_total <- sum(base_counts)
  base_freq <- if (base_total > 0) base_counts / base_total else base_counts

  observed_counts <- count_kmers(analysis_seq, 2L, "^[ACGT]{2}$")
  observed_total <- sum(observed_counts)

  out <- data.frame(
    species = meta$species,
    seq_id = meta$seq_id,
    window_index = meta$window_index,
    start = meta$start,
    end = meta$end,
    window_length = meta$window_length,
    acgt_fraction = meta$acgt_fraction,
    dimer = dimers,
    observed_count = as.integer(observed_counts[dimers]),
    stringsAsFactors = FALSE
  )
  out$observed_count[is.na(out$observed_count)] <- 0L
  out$observed_frequency <- if (observed_total > 0) out$observed_count / observed_total else 0

  split_dimers <- strsplit(out$dimer, "", fixed = TRUE)
  expected_probability <- vapply(
    split_dimers,
    function(x) unname(base_freq[[x[1L]]] * base_freq[[x[2L]]]),
    numeric(1)
  )
  out$expected_probability <- expected_probability
  out$expected_count <- expected_probability * observed_total
  out$ratio_observed_expected <- ifelse(out$expected_count > 0, out$observed_count / out$expected_count, NA_real_)
  out$log2_ratio_observed_expected <- log2_or_na(out$ratio_observed_expected)
  out
}

transition_matrix <- function(seq) {
  bases <- c("A", "C", "G", "T")
  dimers <- count_kmers(seq, 2L, "^[ACGT]{2}$")
  mat <- matrix(0, nrow = 4L, ncol = 4L, dimnames = list(bases, bases))
  if (length(dimers) > 0) {
    for (dimer in names(dimers)) {
      from <- substr(dimer, 1L, 1L)
      to <- substr(dimer, 2L, 2L)
      mat[from, to] <- dimers[[dimer]]
    }
  }
  row_totals <- rowSums(mat)
  probs <- mat
  probs[] <- 0
  valid_rows <- row_totals > 0
  probs[valid_rows, ] <- mat[valid_rows, , drop = FALSE] / row_totals[valid_rows]
  probs
}

markov_hexamer_probability <- function(hexamer, base_freq, transitions) {
  chars <- strsplit(hexamer, "", fixed = TRUE)[[1]]
  prob <- unname(base_freq[[chars[1L]]])
  if (is.na(prob) || prob == 0) {
    return(0)
  }
  for (i in 2:length(chars)) {
    step_prob <- transitions[chars[i - 1L], chars[i]]
    if (is.na(step_prob) || step_prob == 0) {
      return(0)
    }
    prob <- prob * step_prob
  }
  prob
}

bhlh_table_for_window <- function(window_seq, meta, include_revcomp) {
  seq <- toupper(window_seq)
  analysis_seq <- if (include_revcomp) {
    paste0(seq, "NNNNNN", reverse_complement(seq))
  } else {
    seq
  }

  internal_dimers <- as.vector(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0))
  hexamers <- paste0("CA", internal_dimers, "TG")
  bases <- strsplit(analysis_seq, "", fixed = TRUE)[[1]]
  acgt_bases <- bases[bases %in% c("A", "C", "G", "T")]
  base_counts <- table(factor(acgt_bases, levels = c("A", "C", "G", "T")))
  base_total <- sum(base_counts)
  base_freq <- if (base_total > 0) base_counts / base_total else base_counts

  observed_hexamers <- count_kmers(analysis_seq, 6L, "^[ACGT]{6}$")
  observed_valid_hexamer_total <- sum(observed_hexamers)
  observed_bhlh_total <- sum(observed_hexamers[hexamers], na.rm = TRUE)
  transitions <- transition_matrix(analysis_seq)

  expected_probability <- vapply(
    hexamers,
    markov_hexamer_probability,
    numeric(1),
    base_freq = base_freq,
    transitions = transitions
  )

  out <- data.frame(
    species = meta$species,
    seq_id = meta$seq_id,
    window_index = meta$window_index,
    start = meta$start,
    end = meta$end,
    window_length = meta$window_length,
    acgt_fraction = meta$acgt_fraction,
    motif = hexamers,
    internal_dimer = internal_dimers,
    observed_count = as.integer(observed_hexamers[hexamers]),
    stringsAsFactors = FALSE
  )
  out$observed_count[is.na(out$observed_count)] <- 0L
  out$observed_frequency_all_hexamers <- if (observed_valid_hexamer_total > 0) {
    out$observed_count / observed_valid_hexamer_total
  } else {
    0
  }
  out$observed_frequency_bhlh <- if (observed_bhlh_total > 0) {
    out$observed_count / observed_bhlh_total
  } else {
    0
  }
  out$expected_probability_markov <- expected_probability
  out$expected_count_markov <- expected_probability * observed_valid_hexamer_total
  out$ratio_observed_expected <- ifelse(
    out$expected_count_markov > 0,
    out$observed_count / out$expected_count_markov,
    NA_real_
  )
  out$log2_ratio_observed_expected <- log2_or_na(out$ratio_observed_expected)
  out
}

functional_group_from_motif <- function(motif) {
  first_monomer <- substr(motif, 1L, 3L)
  second_monomer <- paste0("CA", chartr("ACGT", "TGCA", substr(motif, 4L, 4L)))
  paste(sort(c(first_monomer, second_monomer), decreasing = TRUE), collapse = "-")
}

functional_bhlh_table_from_sequence_table <- function(sequence_table, meta) {
  sequence_table$functional_group <- vapply(
    sequence_table$motif,
    functional_group_from_motif,
    character(1)
  )

  groups <- sort(unique(sequence_table$functional_group))
  out <- do.call(
    rbind,
    lapply(groups, function(group) {
      rows <- sequence_table[sequence_table$functional_group == group, , drop = FALSE]
      monomers <- strsplit(group, "-", fixed = TRUE)[[1]]
      observed_count <- sum(rows$observed_count)
      expected_count <- sum(rows$expected_count_markov)
      observed_valid_hexamer_total <- unique(ifelse(
        rows$observed_frequency_all_hexamers > 0,
        rows$observed_count / rows$observed_frequency_all_hexamers,
        NA_real_
      ))
      observed_valid_hexamer_total <- observed_valid_hexamer_total[!is.na(observed_valid_hexamer_total)]
      observed_valid_hexamer_total <- if (length(observed_valid_hexamer_total) > 0) {
        observed_valid_hexamer_total[[1]]
      } else {
        0
      }

      data.frame(
        species = meta$species,
        seq_id = meta$seq_id,
        window_index = meta$window_index,
        start = meta$start,
        end = meta$end,
        window_length = meta$window_length,
        acgt_fraction = meta$acgt_fraction,
        functional_group = group,
        monomer_1 = monomers[[1]],
        monomer_2 = monomers[[2]],
        sequence_motifs = paste(rows$motif, collapse = ","),
        n_sequence_motifs = nrow(rows),
        observed_count = observed_count,
        observed_frequency_all_hexamers = if (observed_valid_hexamer_total > 0) {
          observed_count / observed_valid_hexamer_total
        } else {
          0
        },
        expected_probability_markov = sum(rows$expected_probability_markov),
        expected_count_markov = expected_count,
        ratio_observed_expected = if (expected_count > 0) observed_count / expected_count else NA_real_,
        log2_ratio_observed_expected = if (expected_count > 0) log2(observed_count / expected_count) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

analyse_one_window <- function(seq, seq_id, window_row, opts) {
  window_seq <- substr(seq, window_row$start, window_row$end)
  bases <- strsplit(toupper(window_seq), "", fixed = TRUE)[[1]]
  acgt_fraction <- if (length(bases) > 0) {
    mean(bases %in% c("A", "C", "G", "T"))
  } else {
    0
  }

  meta <- list(
    species = opts$species,
    seq_id = seq_id,
    window_index = window_row$window_index,
    start = window_row$start,
    end = window_row$end,
    window_length = nchar(window_seq),
    acgt_fraction = acgt_fraction
  )

  if (acgt_fraction < opts$min_acgt_fraction) {
    return(NULL)
  }

  bhlh <- bhlh_table_for_window(window_seq, meta, opts$include_revcomp)

  list(
    dinucleotides = dinucleotide_table_for_window(window_seq, meta, opts$include_revcomp),
    bhlh = bhlh,
    functional_bhlh = functional_bhlh_table_from_sequence_table(bhlh, meta)
  )
}

append_table <- function(path, table) {
  write.table(
    table,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    append = file.exists(path),
    col.names = !file.exists(path)
  )
}

write_functional_summary_tables <- function(functional_bhlh_path, opts) {
  if (!file.exists(functional_bhlh_path)) {
    return(invisible(NULL))
  }

  functional <- read.delim(functional_bhlh_path, stringsAsFactors = FALSE, check.names = FALSE)
  functional <- functional[!is.na(functional$log2_ratio_observed_expected), , drop = FALSE]
  if (nrow(functional) == 0L) {
    return(invisible(NULL))
  }

  add_direction <- function(table) {
    table$representation <- ifelse(table$log2_ratio_observed_expected >= 0, "over", "under")
    table
  }

  top_n_rows <- do.call(
    rbind,
    lapply(split(functional, functional$functional_group), function(group_table) {
      over <- group_table[group_table$log2_ratio_observed_expected > 0, , drop = FALSE]
      under <- group_table[group_table$log2_ratio_observed_expected < 0, , drop = FALSE]
      over <- over[order(-over$log2_ratio_observed_expected), , drop = FALSE]
      under <- under[order(under$log2_ratio_observed_expected), , drop = FALSE]
      add_direction(rbind(head(over, opts$top_n), head(under, opts$top_n)))
    })
  )
  rownames(top_n_rows) <- NULL
  top_n_path <- file.path(
    opts$output_dir,
    paste0(opts$species, "_functional_CAN_CAN_top_", opts$top_n, "_over_under_by_group.tsv")
  )
  write.table(top_n_rows, top_n_path, sep = "\t", row.names = FALSE, quote = FALSE)

  threshold_rows <- functional[abs(functional$log2_ratio_observed_expected) >= opts$log2_threshold, , drop = FALSE]
  threshold_rows <- add_direction(threshold_rows)
  threshold_rows <- threshold_rows[order(
    threshold_rows$functional_group,
    -abs(threshold_rows$log2_ratio_observed_expected)
  ), , drop = FALSE]
  threshold_path <- file.path(
    opts$output_dir,
    paste0(opts$species, "_functional_CAN_CAN_log2_threshold_", opts$log2_threshold, ".tsv")
  )
  write.table(threshold_rows, threshold_path, sep = "\t", row.names = FALSE, quote = FALSE)

  top_percent_rows <- do.call(
    rbind,
    lapply(split(functional, functional$functional_group), function(group_table) {
      n_keep <- max(1L, ceiling(nrow(group_table) * opts$top_percent / 100))
      over <- group_table[group_table$log2_ratio_observed_expected > 0, , drop = FALSE]
      under <- group_table[group_table$log2_ratio_observed_expected < 0, , drop = FALSE]
      over <- over[order(-over$log2_ratio_observed_expected), , drop = FALSE]
      under <- under[order(under$log2_ratio_observed_expected), , drop = FALSE]
      add_direction(rbind(head(over, n_keep), head(under, n_keep)))
    })
  )
  rownames(top_percent_rows) <- NULL
  top_percent_label <- gsub(".", "p", as.character(opts$top_percent), fixed = TRUE)
  top_percent_path <- file.path(
    opts$output_dir,
    paste0(opts$species, "_functional_CAN_CAN_top_", top_percent_label, "percent_over_under_by_group.tsv")
  )
  write.table(top_percent_rows, top_percent_path, sep = "\t", row.names = FALSE, quote = FALSE)

  message("Wrote: ", top_n_path)
  message("Wrote: ", threshold_path)
  message("Wrote: ", top_percent_path)
  invisible(NULL)
}

process_sequence <- function(seq_id, seq, opts, dinuc_path, bhlh_path, functional_bhlh_path) {
  seq <- gsub("[[:space:]]", "", seq)
  seq_len <- nchar(seq)
  if (seq_len == 0L) {
    message("Skipping empty sequence: ", seq_id)
    return(invisible(NULL))
  }

  windows <- chunk_windows(seq_len, opts$window_size, opts$step_size)
  message("Processing ", seq_id, " (", seq_len, " bp, ", nrow(windows), " windows)")

  worker <- function(i) {
    analyse_one_window(seq, seq_id, windows[i, ], opts)
  }

  results <- if (opts$cores > 1L && .Platform$OS.type != "windows") {
    mclapply(seq_len(nrow(windows)), worker, mc.cores = opts$cores)
  } else {
    lapply(seq_len(nrow(windows)), worker)
  }
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) {
    return(invisible(NULL))
  }

  dinuc <- do.call(rbind, lapply(results, `[[`, "dinucleotides"))
  bhlh <- do.call(rbind, lapply(results, `[[`, "bhlh"))
  functional_bhlh <- do.call(rbind, lapply(results, `[[`, "functional_bhlh"))
  append_table(dinuc_path, dinuc)
  append_table(bhlh_path, bhlh)
  append_table(functional_bhlh_path, functional_bhlh)
  invisible(NULL)
}

process_fasta <- function(opts) {
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  dinuc_path <- file.path(opts$output_dir, paste0(opts$species, "_dinucleotides_by_window.tsv"))
  bhlh_path <- file.path(opts$output_dir, paste0(opts$species, "_bhlh_CANNTG_by_window.tsv"))
  functional_bhlh_path <- file.path(opts$output_dir, paste0(opts$species, "_bhlh_CAN_CAN_functional_by_window.tsv"))
  if (file.exists(dinuc_path)) file.remove(dinuc_path)
  if (file.exists(bhlh_path)) file.remove(bhlh_path)
  if (file.exists(functional_bhlh_path)) file.remove(functional_bhlh_path)

  con <- open_fasta(opts$fasta)
  on.exit(close(con), add = TRUE)

  current_id <- NULL
  current_seq <- list()
  current_seq_n <- 0L

  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    for (line in lines) {
      if (startsWith(line, ">")) {
        if (!is.null(current_id)) {
          process_sequence(current_id, paste0(unlist(current_seq, use.names = FALSE), collapse = ""), opts, dinuc_path, bhlh_path, functional_bhlh_path)
        }
        current_id <- sanitize_seq_name(line)
        current_seq <- list()
        current_seq_n <- 0L
      } else if (!is.null(current_id)) {
        current_seq_n <- current_seq_n + 1L
        current_seq[[current_seq_n]] <- line
      }
    }
  }

  if (!is.null(current_id)) {
    process_sequence(current_id, paste0(unlist(current_seq, use.names = FALSE), collapse = ""), opts, dinuc_path, bhlh_path, functional_bhlh_path)
  }

  message("Wrote: ", dinuc_path)
  message("Wrote: ", bhlh_path)
  message("Wrote: ", functional_bhlh_path)
  write_functional_summary_tables(functional_bhlh_path, opts)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  process_fasta(opts)
}

main()
