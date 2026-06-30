#!/bin/bash
#SBATCH --job-name=functional_window_summaries
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/functional_window_summaries-%j.out
#SBATCH --error=logs/functional_window_summaries-%j.err

set -euo pipefail

module load R/4.1.2

FUNCTIONAL_TABLE=${1:?Usage: sbatch run_functional_summaries.sbatch functional_table.tsv output_dir}
OUTPUT_DIR=${2:?Usage: sbatch run_functional_summaries.sbatch functional_table.tsv output_dir}
TOP_N=${TOP_N:-100}
TOP_PERCENT=${TOP_PERCENT:-1}
LOG2_THRESHOLD=${LOG2_THRESHOLD:-1}
PVALUE_THRESHOLD=${PVALUE_THRESHOLD:-0.05}
MIN_EXPECTED_ZERO=${MIN_EXPECTED_ZERO:-1}
PSEUDOCOUNT=${PSEUDOCOUNT:-0.5}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p logs "$OUTPUT_DIR"

Rscript "$SCRIPT_DIR/summarize_functional_windows.R" \
  --functional_table "$FUNCTIONAL_TABLE" \
  --output_dir "$OUTPUT_DIR" \
  --top_n "$TOP_N" \
  --top_percent "$TOP_PERCENT" \
  --log2_threshold "$LOG2_THRESHOLD" \
  --pvalue_threshold "$PVALUE_THRESHOLD" \
  --min_expected_zero "$MIN_EXPECTED_ZERO" \
  --pseudocount "$PSEUDOCOUNT"
