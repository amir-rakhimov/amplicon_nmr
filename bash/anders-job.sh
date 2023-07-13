#!/usr/bin/env bash
conda activate qiime2
qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filenames-anders.tsv \
  --output-path anders-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data anders-paired-end-demux.qza \
  --o-visualization anders-paired-end-demux.qzv
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences anders-paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences anders-trimmed-paired-end.qza
qiime demux summarize \
  --i-data anders-trimmed-paired-end.qza \
  --o-visualization anders-trimmed-paired-end.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs anders-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 283 \
  --p-trunc-len-r 276 \
  --o-table anders-table-trimmed-dada2-283-276.qza \
  --o-representative-sequences anders-rep-seqs-trimmed-dada2-283-276.qza \
  --o-denoising-stats anders-stats-trimmed-dada2-283-276.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file anders-stats-trimmed-dada2-283-276.qza \
  --o-visualization anders-stats-trimmed-dada2-283-276.qzv
qiime feature-table summarize \
  --i-table anders-table-trimmed-dada2-283-276.qza \
  --o-visualization anders-table-trimmed-dada2-283-276.qzv
qiime feature-table tabulate-seqs \
  --i-data anders-rep-seqs-trimmed-dada2-283-276.qza \
  --o-visualization anders-rep-seqs-trimmed-dada2-283-276.qzv 
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
  --i-reads anders-rep-seqs-trimmed-dada2-283-276.qza \
  --o-classification anders-taxonomy-trimmed-dada2-283-276.qza
qiime metadata tabulate \
  --m-input-file anders-taxonomy-trimmed-dada2-283-276.qza \
  --o-visualization anders-taxonomy-trimmed-dada2-283-276.qzv
