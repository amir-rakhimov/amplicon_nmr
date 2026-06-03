#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
# Load config
source config/bash/config.sh

mkdir -p "${OUTDIR}"
mkdir -p "${qiime2_output_dir}"
mkdir -p "${silva_db_dir}"


job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/qiime2-train-classifier.slurm)
job_id2=$(sbatch --parsable --dependency=afterok:$job_id1 --export=ALL \
            code/slurm-scripts/qiime2-merge-data-and-classify.slurm)
echo "Pipeline submitted"
echo " step 1: ${job_id1}"
echo " step 2: ${job_id2}"
