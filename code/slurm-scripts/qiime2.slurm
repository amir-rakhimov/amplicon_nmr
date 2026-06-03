#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 20
#SBATCH --mem-per-cpu 8g
#SBATCH -J qiime2-pipeline
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
# %x_%j job name and job ID
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
# mkdir -p data/0-raw/fastq/
echo "${RUN_ID}"
bash code/bash-scripts/003-qiime2-pipeline.sh
