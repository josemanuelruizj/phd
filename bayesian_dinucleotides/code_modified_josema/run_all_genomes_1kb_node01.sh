#!/bin/bash

set -euo pipefail

FASTA_DIR=${FASTA_DIR:-/users/genomics/josema/phd/bayesian_dinucleotides/ncbi_species_download/fastas}
OUTDIR=${OUTDIR:-/users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb}
WINDOW_SIZE=${WINDOW_SIZE:-1000}
NODELIST=${NODELIST:-node01}

mkdir -p "$OUTDIR"

for FASTA in \
  "$FASTA_DIR"/*.fna \
  "$FASTA_DIR"/*.fa \
  "$FASTA_DIR"/*.fasta \
  "$FASTA_DIR"/*.fna.gz \
  "$FASTA_DIR"/*.fa.gz \
  "$FASTA_DIR"/*.fasta.gz
do
  [ -e "$FASTA" ] || continue

  sbatch --nodelist="$NODELIST" \
    run_one_genome.sbatch \
    "$FASTA" \
    "$OUTDIR" \
    "$WINDOW_SIZE"
done
