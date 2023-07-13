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
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
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
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 273 \
  --p-trunc-len-r 203 \
  --o-table yasuda-table-trimmed-dada2-273-203.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-273-203.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-273-203.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-273-203.qza \
  --o-visualization yasuda-stats-trimmed-dada2-273-203.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-273-203.qza \
  --o-visualization yasuda-table-trimmed-dada2-273-203.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-273-203.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-273-203.qzv 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 273 \
  --p-trunc-len-r 198 \
  --o-table yasuda-table-trimmed-dada2-273-198.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-273-198.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-273-198.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-273-198.qza \
  --o-visualization yasuda-stats-trimmed-dada2-273-198.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-273-198.qza \
  --o-visualization yasuda-table-trimmed-dada2-273-198.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-273-198.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-273-198.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 265 \
  --p-trunc-len-r 203 \
  --o-table yasuda-table-trimmed-dada2-265-203.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-265-203.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-265-203.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-265-203.qza \
  --o-visualization yasuda-stats-trimmed-dada2-265-203.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-265-203.qza \
  --o-visualization yasuda-table-trimmed-dada2-265-203.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-265-203.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-265-203.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 265 \
  --p-trunc-len-r 198 \
  --o-table yasuda-table-trimmed-dada2-265-198.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-265-198.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-265-198.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-265-198.qza \
  --o-visualization yasuda-stats-trimmed-dada2-265-198.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-265-198.qza \
  --o-visualization yasuda-table-trimmed-dada2-265-198.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-265-198.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-265-198.qzv
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
  --o-reads silva-138-ssu-nr99-seqs-cleaned-341f-806r.qza
qiime rescript dereplicate \
  --i-sequences silva-138-ssu-nr99-seqs-cleaned-341f-806r.qza \
  --i-taxa silva-138-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-ssu-nr99-seqs-cleaned-341f-806r-uniq.qza \
  --o-dereplicated-taxa  silva-138-ssu-nr99-tax-341f-806r-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-ssu-nr99-seqs-cleaned-341f-806r-uniq.qza \
  --i-reference-taxonomy silva-138-ssu-nr99-tax-341f-806r-derep-uniq.qza \
  --o-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-313-229.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-313-229.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-313-229.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-313-229.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-284-229.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-284-229.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-284-229.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-284-229.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-284-203.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-284-203.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-284-203.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-284-203.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-273-203.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-273-203.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-273-203.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-273-203.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-273-198.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-273-198.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-273-198.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-273-198.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-265-203.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-265-203.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-265-203.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-265-203.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-classifier-341f-806r-derep-uniq.qza \
  --i-reads yasuda-rep-seqs-trimmed-dada2-265-198.qza \
  --o-classification yasuda-taxonomy-trimmed-dada2-265-198.qza
qiime metadata tabulate \
  --m-input-file yasuda-taxonomy-trimmed-dada2-265-198.qza \
  --o-visualization yasuda-taxonomy-trimmed-dada2-265-198.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-313-229.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-313-229.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-313-229.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-284-229.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-284-229.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-284-229.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-284-203.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-284-203.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-284-203.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-273-203.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-273-203.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-273-203.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-273-198.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-273-198.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-273-198.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-265-203.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-265-203.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-265-203.qzv
qiime taxa barplot \
  --i-table yasuda-table-trimmed-dada2-265-198.qza \
  --i-taxonomy yasuda-taxonomy-trimmed-dada2-265-198.qza \
  --m-metadata-file filenames-yasuda.tsv \
  --o-visualization yasuda-taxa-bar-plots-trimmed-dada2-265-198.qzv

