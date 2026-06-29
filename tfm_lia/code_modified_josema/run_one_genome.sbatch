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

FASTA=${1:?Usage: sbatch run_one_genome.sbatch genome.fa.gz output_dir [window_size] [overlap]}
OUTPUT_DIR=${2:?Usage: sbatch run_one_genome.sbatch genome.fa.gz output_dir [window_size] [overlap]}
WINDOW_SIZE=${3:-50000}
OVERLAP=${OVERLAP:-${4:-4}}
SPECIES=${SPECIES:-$(basename "$FASTA")}
SPECIES=${SPECIES%.gz}
SPECIES=${SPECIES%.fasta}
SPECIES=${SPECIES%.fa}
SPECIES=${SPECIES%.fna}
SPECIES=${SPECIES%.fas}
STEP_SIZE=${STEP_SIZE:-$((WINDOW_SIZE - OVERLAP))}
TOP_N=${TOP_N:-100}
LOG2_THRESHOLD=${LOG2_THRESHOLD:-1}
TOP_PERCENT=${TOP_PERCENT:-1}
CORES=${SLURM_CPUS_PER_TASK:-1}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$STEP_SIZE" -lt 1 ]; then
  echo "ERROR: step size must be at least 1. Use an overlap smaller than window_size." >&2
  exit 1
fi

mkdir -p logs "$OUTPUT_DIR"

Rscript "$SCRIPT_DIR/genome_window_analysis.R" \
  --fasta "$FASTA" \
  --species "$SPECIES" \
  --output_dir "$OUTPUT_DIR" \
  --window_size "$WINDOW_SIZE" \
  --step_size "$STEP_SIZE" \
  --cores "$CORES" \
  --top_n "$TOP_N" \
  --log2_threshold "$LOG2_THRESHOLD" \
  --top_percent "$TOP_PERCENT" \
  --include_revcomp TRUE
