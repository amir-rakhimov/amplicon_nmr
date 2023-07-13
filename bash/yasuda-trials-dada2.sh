#!/usr/bin/env bash
conda activate qiime2
qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filenames-yasuda.tsv \
  --output-path yasuda-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data yasuda-paired-end-demux.qza \
  --o-visualization yasuda-paired-end-demux.qzv
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences yasuda-paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
  --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences yasuda-trimmed-paired-end.qza
qiime demux summarize \
  --i-data yasuda-trimmed-paired-end.qza \
  --o-visualization yasuda-trimmed-paired-end.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 313 \
  --p-trunc-len-r 229 \
  --o-table yasuda-table-trimmed-dada2-313-229.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-313-229.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-313-229.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-313-229.qza \
  --o-visualization yasuda-stats-trimmed-dada2-313-229.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-313-229.qza \
  --o-visualization yasuda-table-trimmed-dada2-313-229.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-313-229.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-313-229.qzv 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 229 \
  --o-table yasuda-table-trimmed-dada2-284-229.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-284-229.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-284-229.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-284-229.qza \
  --o-visualization yasuda-stats-trimmed-dada2-284-229.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-284-229.qza \
  --o-visualization yasuda-table-trimmed-dada2-284-229.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-284-229.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-284-229.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 203 \
  --o-table yasuda-table-trimmed-dada2-284-203.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-284-203.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-284-203.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-284-203.qza \
  --o-visualization yasuda-stats-trimmed-dada2-284-203.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-284-203.qza \
  --o-visualization yasuda-table-trimmed-dada2-284-203.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-284-203.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-284-203.qzv