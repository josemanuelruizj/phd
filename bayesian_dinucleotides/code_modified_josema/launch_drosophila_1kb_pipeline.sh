#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FASTA=${FASTA:-/users/genomics/josema/phd/bayesian_dinucleotides/ncbi_species_download/fastas/Drosophila_melanogaster__GCF_000001215.4__GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna}
ANNOTATION=${ANNOTATION:-/users/genomics/josema/phd/bayesian_dinucleotides/ncbi_species_download/annotations/Drosophila_melanogaster__GCF_000001215.4__genomic.gtf}
ANNOTATION_DIR=${ANNOTATION_DIR:-$(dirname "$ANNOTATION")}
OUTDIR=${OUTDIR:-/users/genomics/josema/phd/bayesian_dinucleotides/results_drosophila_1kb}

WINDOW_SIZE=${WINDOW_SIZE:-1000}
TOP_N=${TOP_N:-1000}
TOP_PERCENT=${TOP_PERCENT:-3}
LOG2_THRESHOLD=${LOG2_THRESHOLD:-1}
PVALUE_THRESHOLD=${PVALUE_THRESHOLD:-0.1}
SPECIES=${SPECIES:-Drosophila_melanogaster__GCF_000001215.4__GCF_000001215.4_Release_6_plus_ISO1_MT_genomic}

FUNCTIONAL_TABLE="$OUTDIR/${SPECIES}_bhlh_CAN_CAN_functional_by_window.tsv"
SUMMARY_DIR="$OUTDIR/summaries"
CHIPSEEKER_DIR="$OUTDIR/chipseeker_summaries"

mkdir -p "$OUTDIR" "$SUMMARY_DIR" "$CHIPSEEKER_DIR" "$SCRIPT_DIR/logs"

WINDOW_JOB=$(sbatch --parsable \
  "$SCRIPT_DIR/run_one_genome.sbatch" \
  "$FASTA" \
  "$OUTDIR" \
  "$WINDOW_SIZE")

echo "Submitted window analysis job: $WINDOW_JOB"

SUMMARY_JOB=$(TOP_N="$TOP_N" \
  TOP_PERCENT="$TOP_PERCENT" \
  LOG2_THRESHOLD="$LOG2_THRESHOLD" \
  PVALUE_THRESHOLD="$PVALUE_THRESHOLD" \
  sbatch --parsable \
    --dependency="afterok:$WINDOW_JOB" \
    "$SCRIPT_DIR/run_functional_summaries.sbatch" \
    "$FUNCTIONAL_TABLE" \
    "$SUMMARY_DIR")

echo "Submitted summary job: $SUMMARY_JOB after $WINDOW_JOB"

CHIPSEEKER_JOB=$(ANNOTATION_FILE="$ANNOTATION" \
  sbatch --parsable \
    --dependency="afterok:$SUMMARY_JOB" \
    "$SCRIPT_DIR/run_chipseeker_summaries.sbatch" \
    "$SUMMARY_DIR" \
    "$ANNOTATION_DIR" \
    "$CHIPSEEKER_DIR" \
    "$SPECIES")

echo "Submitted ChIPseeker job: $CHIPSEEKER_JOB after $SUMMARY_JOB"
echo ""
echo "Output directory: $OUTDIR"
echo "Summaries: $SUMMARY_DIR"
echo "ChIPseeker summaries: $CHIPSEEKER_DIR"
