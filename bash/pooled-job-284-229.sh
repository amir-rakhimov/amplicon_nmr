#!/usr/bin/env bash
conda activate qiime2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs pooled-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 229 \
  --o-table pooled-table-trimmed-dada2-284-229.qza \
  --o-representative-sequences pooled-rep-seqs-trimmed-dada2-284-229.qza \
  --o-denoising-stats pooled-stats-trimmed-dada2-284-229.qza
qiime metadata tabulate \
  --m-input-file pooled-stats-trimmed-dada2-284-229.qza \
  --o-visualization pooled-stats-trimmed-dada2-284-229.qzv
qiime feature-table summarize \
  --i-table pooled-table-trimmed-dada2-284-229.qza \
  --o-visualization pooled-table-trimmed-dada2-284-229.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-284-229.qza \
  --o-visualization pooled-rep-seqs-trimmed-dada2-284-229.qzv
qiime feature-table filter-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-284-229.qza \
  --m-metadata-file pooled-rep-seqs-trimmed-dada2-284-229.qza \
  --p-where 'length(sequence) > 299' \
  --o-filtered-data pooled-filtered-rep-seqs-trimmed-dada2-284-229.qza
qiime feature-table filter-features \
 --i-table pooled-table-trimmed-dada2-284-229.qza \
 --m-metadata-file pooled-filtered-rep-seqs-trimmed-dada2-284-229.qza \
 --p-no-exclude-ids \
 --o-filtered-table pooled-filtered-table-trimmed-dada2-284-229.qza
qiime feature-table summarize \
  --i-table pooled-filtered-table-trimmed-dada2-284-229.qza \
  --o-visualization pooled-filtered-table-trimmed-dada2-284-229.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-filtered-rep-seqs-trimmed-dada2-284-229.qza \
  --o-visualization pooled-filtered-rep-seqs-trimmed-dada2-284-229.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-1-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads pooled-filtered-rep-seqs-trimmed-dada2-284-229.qza \
  --o-classification pooled-taxonomy-trimmed-dada2-284-229.qza
qiime metadata tabulate \
  --m-input-file pooled-taxonomy-trimmed-dada2-284-229.qza \
  --o-visualization pooled-taxonomy-trimmed-dada2-284-229.qzv
qiime taxa barplot \
  --i-table pooled-filtered-table-trimmed-dada2-284-229.qza \
  --i-taxonomy pooled-taxonomy-trimmed-dada2-284-229.qza \
  --m-metadata-file filenames-pooled-raw-supercomp.tsv \
  --o-visualization pooled-taxa-bar-plots-trimmed-dada2-284-229.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences pooled-filtered-rep-seqs-trimmed-dada2-284-229.qza \
  --o-alignment pooled-aligned-rep-seqs-trimmed-dada2-284-229.qza \
  --o-masked-alignment pooled-masked-aligned-rep-seqs-trimmed-dada2-284-229.qza \
  --o-tree pooled-unrooted-tree-trimmed-dada2-284-229.qza \
  --o-rooted-tree pooled-rooted-tree-trimmed-dada2-284-229.qza