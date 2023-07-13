#!/usr/bin/env bash
conda activate qiime2
qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filenames-pooled-raw-supercomp.tsv \
  --output-path pooled-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data pooled-paired-end-demux.qza \
  --o-visualization pooled-paired-end-demux.qzv
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences pooled-paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences pooled-trimmed-paired-end.qza
qiime demux summarize \
  --i-data pooled-trimmed-paired-end.qza \
  --o-visualization pooled-trimmed-paired-end.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs pooled-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 203 \
  --o-table pooled-table-trimmed-dada2-284-203.qza \
  --o-representative-sequences pooled-rep-seqs-trimmed-dada2-284-203.qza \
  --o-denoising-stats pooled-stats-trimmed-dada2-284-203.qza
qiime metadata tabulate \
  --m-input-file pooled-stats-trimmed-dada2-284-203.qza \
  --o-visualization pooled-stats-trimmed-dada2-284-203.qzv
qiime feature-table summarize \
  --i-table pooled-table-trimmed-dada2-284-203.qza \
  --o-visualization pooled-table-trimmed-dada2-284-203.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-284-203.qza \
  --o-visualization pooled-rep-seqs-trimmed-dada2-284-203.qzv
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138-1-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138-1-ssu-nr99-tax.qza
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138-1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138-1-ssu-nr99-seqs.qza
qiime rescript cull-seqs \
  --i-sequences silva-138-1-ssu-nr99-seqs.qza \
  --o-clean-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza
qiime feature-classifier extract-reads \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-138-1-ssu-nr99-seqs-cleaned-341f-805r.qza
qiime rescript dereplicate \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned-341f-805r.qza \
  --i-taxa silva-138-1-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-1-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --o-dereplicated-taxa  silva-138-1-ssu-nr99-tax-341f-805r-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-1-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --i-reference-taxonomy silva-138-1-ssu-nr99-tax-341f-805r-derep-uniq.qza \
  --o-classifier silva-138-1-ssu-nr99-classifier-341f-805r-derep-uniq.qza
qiime feature-table filter-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-284-203.qza \
  --m-metadata-file pooled-rep-seqs-trimmed-dada2-284-203.qza \
  --p-where 'length(sequence) > 299' \
  --o-filtered-data pooled-filtered-rep-seqs-trimmed-dada2-284-203.qza
qiime feature-table filter-features \
 --i-table pooled-table-trimmed-dada2-284-203.qza \
 --m-metadata-file pooled-filtered-rep-seqs-trimmed-dada2-284-203.qza \
 --p-no-exclude-ids \
 --o-filtered-table pooled-filtered-table-trimmed-dada2-284-203.qza
qiime feature-table summarize \
  --i-table pooled-filtered-table-trimmed-dada2-284-203.qza \
  --o-visualization pooled-filtered-table-trimmed-dada2-284-203.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-filtered-rep-seqs-trimmed-dada2-284-203.qza \
  --o-visualization pooled-filtered-rep-seqs-trimmed-dada2-284-203.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-1-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads pooled-filtered-rep-seqs-trimmed-dada2-284-203.qza \
  --o-classification pooled-taxonomy-trimmed-dada2-284-203.qza
qiime metadata tabulate \
  --m-input-file pooled-taxonomy-trimmed-dada2-284-203.qza \
  --o-visualization pooled-taxonomy-trimmed-dada2-284-203.qzv
qiime taxa barplot \
  --i-table pooled-filtered-table-trimmed-dada2-284-203.qza \
  --i-taxonomy pooled-taxonomy-trimmed-dada2-284-203.qza \
  --m-metadata-file filenames-pooled-raw-supercomp.tsv \
  --o-visualization pooled-taxa-bar-plots-trimmed-dada2-284-203.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences pooled-filtered-rep-seqs-trimmed-dada2-284-203.qza \
  --o-alignment pooled-aligned-rep-seqs-trimmed-dada2-284-203.qza \
  --o-masked-alignment pooled-masked-aligned-rep-seqs-trimmed-dada2-284-203.qza \
  --o-tree pooled-unrooted-tree-trimmed-dada2-284-203.qza \
  --o-rooted-tree pooled-rooted-tree-trimmed-dada2-284-203.qza