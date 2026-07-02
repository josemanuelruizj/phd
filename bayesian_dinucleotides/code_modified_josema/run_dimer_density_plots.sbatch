#!/bin/bash
#SBATCH --job-name=dimer_density_plots
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dimer_density_plots-%j.out
#SBATCH --error=logs/dimer_density_plots-%j.err

set -euo pipefail

module load R/4.3.2

INPUT_DIR=${1:-/users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb}
METHYLATION_TABLE=${2:-/users/genomics/josema/phd/bayesian_dinucleotides/code_modified_josema/dnmt_methylation_info_download_species_UPDATED.tsv}
OUTPUT_DIR=${3:-$INPUT_DIR/dimer_density_plots}
X_MIN=${X_MIN:-0}
X_MAX=${X_MAX:-2.2}
GRID_N=${GRID_N:-512}
RATIO_COL=${RATIO_COL:-ratio_observed_expected}
BW=${BW:-nrd0}
SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
SCRIPT_PATH="$SCRIPT_DIR/plot_dimer_ratio_densities.R"

if [ ! -f "$SCRIPT_PATH" ]; then
  echo "ERROR: Cannot find plot_dimer_ratio_densities.R in $SCRIPT_DIR" >&2
  echo "Run sbatch from the code_modified_josema directory or set SLURM_SUBMIT_DIR correctly." >&2
  exit 1
fi

mkdir -p logs "$OUTPUT_DIR"

Rscript "$SCRIPT_PATH" \
  --input_dir "$INPUT_DIR" \
  --methylation_table "$METHYLATION_TABLE" \
  --output_dir "$OUTPUT_DIR" \
  --ratio_col "$RATIO_COL" \
  --x_min "$X_MIN" \
  --x_max "$X_MAX" \
  --grid_n "$GRID_N" \
  --bw "$BW"
