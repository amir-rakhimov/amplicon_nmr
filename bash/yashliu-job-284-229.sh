#!/usr/bin/env bash
conda activate qiime2
qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filenames-yashliu-raw-supercomp.tsv \
  --output-path yashliu-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data yashliu-paired-end-demux.qza \
  --o-visualization yashliu-paired-end-demux.qzv
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences yashliu-paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences yashliu-trimmed-paired-end.qza
qiime demux summarize \
  --i-data yashliu-trimmed-paired-end.qza \
  --o-visualization yashliu-trimmed-paired-end.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yashliu-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 229 \
  --o-table yashliu-table-trimmed-dada2-284-229.qza \
  --o-representative-sequences yashliu-rep-seqs-trimmed-dada2-284-229.qza \
  --o-denoising-stats yashliu-stats-trimmed-dada2-284-229.qza
qiime metadata tabulate \
  --m-input-file yashliu-stats-trimmed-dada2-284-229.qza \
  --o-visualization yashliu-stats-trimmed-dada2-284-229.qzv
qiime feature-table summarize \
  --i-table yashliu-table-trimmed-dada2-284-229.qza \
  --o-visualization yashliu-table-trimmed-dada2-284-229.qzv
qiime feature-table tabulate-seqs \
  --i-data yashliu-rep-seqs-trimmed-dada2-284-229.qza \
  --o-visualization yashliu-rep-seqs-trimmed-dada2-284-229.qzv
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
qiime rescript filter-seqs-length \
  --i-sequences yashliu-rep-seqs-trimmed-dada2-284-229.qza \
  --p-global-min 300 \
  --o-filtered-seqs yashliu-filtered-rep-seqs-trimmed-dada2-284-229.qza \
  --o-discarded-seqs yashliu-discarded-seqs-trimmed-dada2-284-229.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads yashliu-filtered-rep-seqs-trimmed-dada2-284-229.qza \
  --o-classification yashliu-taxonomy-filtered-trimmed-dada2-284-229.qza
qiime metadata tabulate \
  --m-input-file yashliu-taxonomy-filtered-trimmed-dada2-284-229.qza \
  --o-visualization yashliu-taxonomy-filtered-trimmed-dada2-284-229.qzv
qiime taxa barplot \
  --i-table yashliu-table-trimmed-dada2-284-229.qza \
  --i-taxonomy yashliu-taxonomy-filtered-trimmed-dada2-284-229.qza \
  --m-metadata-file filenames-yashliu-raw-supercomp.tsv \
  --o-visualization yashliu-taxa-bar-plots-filtered-trimmed-dada2-284-229.qzv