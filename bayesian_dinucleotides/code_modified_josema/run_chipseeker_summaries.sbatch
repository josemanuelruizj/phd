#!/bin/bash
#SBATCH --job-name=functional_chipseeker
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=logs/functional_chipseeker-%j.out
#SBATCH --error=logs/functional_chipseeker-%j.err

set -euo pipefail

module load R/4.3.2

SUMMARY_DIR=${1:?Usage: sbatch run_chipseeker_summaries.sbatch summary_dir annotation_dir output_dir [species]}
ANNOTATION_DIR=${2:?Usage: sbatch run_chipseeker_summaries.sbatch summary_dir annotation_dir output_dir [species]}
OUTPUT_DIR=${3:?Usage: sbatch run_chipseeker_summaries.sbatch summary_dir annotation_dir output_dir [species]}
SPECIES=${4:-${SPECIES:-}}
FASTA_DIR=${FASTA_DIR:-}
ANNOTATION_FILE=${ANNOTATION_FILE:-}
TSS_UPSTREAM=${TSS_UPSTREAM:-3000}
TSS_DOWNSTREAM=${TSS_DOWNSTREAM:-3000}
SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
SCRIPT_PATH="$SCRIPT_DIR/annotate_functional_summaries_chipseeker.R"

if [ ! -f "$SCRIPT_PATH" ]; then
  echo "ERROR: Cannot find annotate_functional_summaries_chipseeker.R in $SCRIPT_DIR" >&2
  echo "Run sbatch from the code_modified_josema directory or set SLURM_SUBMIT_DIR correctly." >&2
  exit 1
fi

mkdir -p logs "$OUTPUT_DIR"

ARGS=(
  --summary_dir "$SUMMARY_DIR"
  --annotation_dir "$ANNOTATION_DIR"
  --output_dir "$OUTPUT_DIR"
  --tss_upstream "$TSS_UPSTREAM"
  --tss_downstream "$TSS_DOWNSTREAM"
)

if [ -n "$SPECIES" ]; then
  ARGS+=(--species "$SPECIES")
fi

if [ -n "$FASTA_DIR" ]; then
  ARGS+=(--fasta_dir "$FASTA_DIR")
fi

if [ -n "$ANNOTATION_FILE" ]; then
  ARGS+=(--annotation_file "$ANNOTATION_FILE")
fi

Rscript "$SCRIPT_PATH" "${ARGS[@]}"
