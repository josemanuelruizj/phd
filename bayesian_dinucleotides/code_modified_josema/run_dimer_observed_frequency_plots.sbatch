#!/bin/bash
#SBATCH --job-name=dimer_freq_plots
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G
#SBATCH --time=06:00:00
#SBATCH --output=logs/dimer_freq_plots-%j.out
#SBATCH --error=logs/dimer_freq_plots-%j.err

set -euo pipefail

module load R/4.3.2

INPUT_DIR=${1:-/users/genomics/josema/phd/bayesian_dinucleotides/results_all_genomes_1kb}
METHYLATION_TABLE=${2:-/users/genomics/josema/phd/bayesian_dinucleotides/code_modified_josema/dnmt_methylation_info_download_species_UPDATED.tsv}
OUTPUT_DIR=${3:-$INPUT_DIR/dimer_observed_frequency_plots}
FREQUENCY_COL=${FREQUENCY_COL:-observed_frequency}
X_MAX=${X_MAX:-NA}
SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
SCRIPT_PATH="$SCRIPT_DIR/plot_dimer_observed_frequency.R"

if [ ! -f "$SCRIPT_PATH" ]; then
  echo "ERROR: Cannot find plot_dimer_observed_frequency.R in $SCRIPT_DIR" >&2
  echo "Run sbatch from the code_modified_josema directory or set SLURM_SUBMIT_DIR correctly." >&2
  exit 1
fi

mkdir -p logs "$OUTPUT_DIR"

Rscript "$SCRIPT_PATH" \
  --input_dir "$INPUT_DIR" \
  --methylation_table "$METHYLATION_TABLE" \
  --output_dir "$OUTPUT_DIR" \
  --frequency_col "$FREQUENCY_COL" \
  --x_max "$X_MAX"
