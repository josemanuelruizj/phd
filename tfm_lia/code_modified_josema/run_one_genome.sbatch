#!/bin/bash
#SBATCH --job-name=genome_window_analysis
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --output=logs/genome_window_analysis-%j.out
#SBATCH --error=logs/genome_window_analysis-%j.err

set -euo pipefail

module load R/4.1.2

FASTA=${1:?Usage: sbatch run_one_genome.sbatch genome.fa.gz species output_dir [window_size]}
SPECIES=${2:?Usage: sbatch run_one_genome.sbatch genome.fa.gz species output_dir [window_size]}
OUTPUT_DIR=${3:?Usage: sbatch run_one_genome.sbatch genome.fa.gz species output_dir [window_size]}
WINDOW_SIZE=${4:-50000}
STEP_SIZE=${STEP_SIZE:-$WINDOW_SIZE}
CORES=${SLURM_CPUS_PER_TASK:-1}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p logs "$OUTPUT_DIR"

Rscript "$SCRIPT_DIR/genome_window_analysis.R" \
  --fasta "$FASTA" \
  --species "$SPECIES" \
  --output_dir "$OUTPUT_DIR" \
  --window_size "$WINDOW_SIZE" \
  --step_size "$STEP_SIZE" \
  --cores "$CORES" \
  --include_revcomp TRUE
