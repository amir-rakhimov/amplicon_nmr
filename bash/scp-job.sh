#!/usr/bin/env bash
-S /bin/bash -V -cwd -l short -l s_vmem=8G, mem_req=8G -l d_rt=24:00:00, s_rt=24:00:00 -pe def_slot 4


conda init bash
conda activate /home/rakhimov/miniconda3/envs/qiime2-2022.8 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs song-paired-end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table song-table-dada2-no-trunc.qza \
  --o-representative-sequences song-rep-seqs-dada2-no-trunc.qza \
  --o-denoising-stats song-stats-dada2-no-trunc.qza \
  --p-n-threads 0