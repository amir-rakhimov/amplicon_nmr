#!/usr/bin/env bash
conda activate qiime2
qiime tools import \
  --type  'SampleData[SequencesWithQuality]' \
  --input-path filenames-biagi.tsv \
  --output-path biagi-single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data biagi-single-end-demux.qza \
  --o-visualization biagi-single-end-demux.qzv
qiime cutadapt trim-single \
  --i-demultiplexed-sequences biagi-single-end-demux.qza \
  --p-cores 8 \
  --p-front CCTACGGGNGGCWGCAG \
  --p-adapter GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences biagi-trimmed-single-end.qza
qiime demux summarize \
  --i-data biagi-trimmed-single-end.qza \
  --o-visualization biagi-trimmed-single-end.qzv
qiime dada2 denoise-single \
  --i-demultiplexed-seqs biagi-trimmed-single-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences biagi-rep-seqs-trimmed-dada2-no-trunc.qza \
  --o-table biagi-table-trimmed-dada2-no-trunc.qza \
  --o-denoising-stats biagi-stats-trimmed-dada2-no-trunc.qza \
  --p-n-threads 0  
qiime metadata tabulate \
  --m-input-file biagi-stats-trimmed-dada2-no-trunc.qza \
  --o-visualization biagi-stats-trimmed-dada2-no-trunc.qzv
qiime feature-table summarize \
  --i-table biagi-table-trimmed-dada2-no-trunc.qza \
  --o-visualization biagi-table-trimmed-dada2-no-trunc.qzv
qiime feature-table tabulate-seqs \
  --i-data biagi-rep-seqs-trimmed-dada2-no-trunc.qza \
  --o-visualization biagi-rep-seqs-trimmed-dada2-no-trunc.qzv 
qiime rescript get-silva-data \
    --p-version '138' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138-ssu-nr99-tax.qza
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138-ssu-nr99-seqs.qza
qiime rescript cull-seqs \
  --i-sequences silva-138-ssu-nr99-seqs.qza \
  --o-clean-sequences silva-138-ssu-nr99-seqs-cleaned.qza
qiime feature-classifier extract-reads \
  --i-sequences silva-138-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-138-ssu-nr99-seqs-cleaned-341f-805r.qza
qiime rescript dereplicate \
  --i-sequences silva-138-ssu-nr99-seqs-cleaned-341f-805r.qza \
  --i-taxa silva-138-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --o-dereplicated-taxa  silva-138-ssu-nr99-tax-341f-805r-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --i-reference-taxonomy silva-138-ssu-nr99-tax-341f-805r-derep-uniq.qza \
  --o-classifier silva-138-ssu-nr99-classifier-341f-805r-derep-uniq.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads biagi-rep-seqs-trimmed-dada2-no-trunc.qza \
  --o-classification biagi-taxonomy-trimmed-dada2-no-trunc.qza
qiime metadata tabulate \
  --m-input-file biagi-taxonomy-trimmed-dada2-no-trunc.qza \
  --o-visualization biagi-taxonomy-trimmed-dada2-no-trunc.qzv
qiime taxa barplot \
  --i-table biagi-table-trimmed-dada2-no-trunc.qza \
  --i-taxonomy biagi-taxonomy-trimmed-dada2-no-trunc.qza \
  --m-metadata-file filenames-biagi.tsv \
  --o-visualization biagi-taxa-bar-plots-trimmed-dada2-no-trunc.qzv