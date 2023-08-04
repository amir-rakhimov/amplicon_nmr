#!/usr/bin/env bash
trunc_len_f=284
trunc_len_r=203
front_primer=CCTACGGGNGGCWGCAG
rev_primer=GACTACHVGGGTATCTAATCC
front_primer_position=341f
rev_primer_position=805r
ref_db_trunc_len=465
metadata_file="filenames-pooled-raw-supercomp.tsv"
conda activate qiime2
qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${metadata_file} \
  --output-path pooled-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data pooled-paired-end-demux.qza \
  --o-visualization pooled-paired-end-demux.qzv
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences pooled-paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f ${front_primer} \
  --p-front-r ${rev_primer} \
  --o-trimmed-sequences pooled-trimmed-paired-end.qza
qiime demux summarize \
  --i-data pooled-trimmed-paired-end.qza \
  --o-visualization pooled-trimmed-paired-end.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs pooled-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f ${trunc_len_f} \
  --p-trunc-len-r ${trunc_len_r} \
  --o-table pooled-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-representative-sequences pooled-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-denoising-stats pooled-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime metadata tabulate \
  --m-input-file pooled-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization pooled-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-table summarize \
  --i-table pooled-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization pooled-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
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
  --p-f-primer ${front_primer} \
  --p-r-primer ${rev_primer} \
  --p-trunc-len ${ref_db_trunc_len} \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}.qza
qiime rescript dereplicate \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}.qza \
  --i-taxa silva-138-1-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}-uniq.qza \
  --o-dereplicated-taxa  silva-138-1-ssu-nr99-tax-${front_primer_position}-${rev_primer_position}-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}-uniq.qza \
  --i-reference-taxonomy silva-138-1-ssu-nr99-tax-${front_primer_position}-${rev_primer_position}-derep-uniq.qza \
  --o-classifier silva-138-1-ssu-nr99-classifier-${front_primer_position}-${rev_primer_position}-derep-uniq.qza
qiime feature-table filter-seqs \
  --i-data pooled-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-metadata-file pooled-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-where 'length(sequence) > 299' \
  --o-filtered-data pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime feature-table filter-features \
 --i-table pooled-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
 --m-metadata-file pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
 --p-no-exclude-ids \
 --o-filtered-table pooled-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime feature-table summarize \
  --i-table pooled-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-1-ssu-nr99-classifier-${front_primer_position}-${rev_primer_position}-derep-uniq.qza \
  --i-reads pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-classification pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime taxa filter-table \
  --i-table pooled-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --i-taxonomy pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-table pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza
qiime feature-table summarize \
  --i-table pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qzv
mv pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza pooled-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime taxa filter-seqs \
  --i-sequences pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --i-taxonomy pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-sequences pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza
qiime feature-table tabulate-seqs \
  --i-data pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza \
  --o-visualization pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qzv
mv pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime metadata tabulate \
  --m-input-file pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime taxa barplot \
  --i-table pooled-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --i-taxonomy pooled-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-metadata-file ${metadata_file} \
  --o-visualization pooled-taxa-bar-plots-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences pooled-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-alignment pooled-aligned-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-masked-alignment pooled-masked-aligned-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-tree pooled-unrooted-tree-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-rooted-tree pooled-rooted-tree-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza