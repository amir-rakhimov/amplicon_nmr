#!/usr/bin/env bash
conda activate qiime2
qiime tools import \
  --type  'SampleData[SequencesWithQuality]' \
  --input-path filenames-single-pooled-raw-supercomp.tsv \
  --output-path pooled-single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data pooled-single-end-demux.qza \
  --o-visualization pooled-single-end-demux.qzv
qiime cutadapt trim-single \
  --i-demultiplexed-sequences pooled-single-end-demux.qza \
  --p-cores 8 \
  --p-front CCTACGGGNGGCWGCAG \
  --o-trimmed-sequences pooled-trimmed-single-end.qza
qiime demux summarize \
  --i-data pooled-trimmed-single-end.qza \
  --o-visualization pooled-trimmed-single-end.qzv
qiime dada2 denoise-single \
  --i-demultiplexed-seqs pooled-trimmed-single-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 234 \
  --o-table pooled-table-trimmed-dada2-234.qza \
  --o-representative-sequences pooled-rep-seqs-trimmed-dada2-234.qza \
  --o-denoising-stats pooled-stats-trimmed-dada2-234.qza
qiime metadata tabulate \
  --m-input-file pooled-stats-trimmed-dada2-234.qza \
  --o-visualization pooled-stats-trimmed-dada2-234.qzv
qiime feature-table summarize \
  --i-table pooled-table-trimmed-dada2-234.qza \
  --o-visualization pooled-table-trimmed-dada2-234.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-234.qza \
  --o-visualization pooled-rep-seqs-trimmed-dada2-234.qzv
qiime feature-classifier extract-reads \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 234 \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-138-1-ssu-nr99-seqs-cleaned-341f-234.qza
qiime rescript dereplicate \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned-341f-234.qza \
  --i-taxa silva-138-1-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-1-ssu-nr99-seqs-cleaned-341f-234-uniq.qza \
  --o-dereplicated-taxa  silva-138-1-ssu-nr99-tax-341f-234-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-1-ssu-nr99-seqs-cleaned-341f-234-uniq.qza \
  --i-reference-taxonomy silva-138-1-ssu-nr99-tax-341f-234-derep-uniq.qza \
  --o-classifier silva-138-1-ssu-nr99-classifier-341f-234-derep-uniq.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-1-ssu-nr99-classifier-341f-234-derep-uniq.qza \
  --i-reads pooled-rep-seqs-trimmed-dada2-234.qza \
  --o-classification pooled-taxonomy-trimmed-dada2-234.qza
qiime metadata tabulate \
  --m-input-file pooled-taxonomy-trimmed-dada2-234.qza \
  --o-visualization pooled-taxonomy-trimmed-dada2-234.qzv
qiime taxa barplot \
  --i-table pooled-table-trimmed-dada2-234.qza \
  --i-taxonomy pooled-taxonomy-trimmed-dada2-234.qza \
  --m-metadata-file filenames-single-pooled-raw-supercomp.tsv \
  --o-visualization pooled-taxa-bar-plots-trimmed-dada2-234.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences pooled-rep-seqs-trimmed-dada2-234.qza \
  --o-alignment pooled-aligned-rep-seqs-trimmed-dada2-234.qza \
  --o-masked-alignment pooled-masked-aligned-rep-seqs-trimmed-dada2-234.qza \
  --o-tree pooled-unrooted-tree-trimmed-dada2-234.qza \
  --o-rooted-tree pooled-rooted-tree-trimmed-dada2-234.qza