#!/usr/bin/env bash
author_name=biagi
trunc_len_f=0
trunc_len_r=0
front_primer=CCTACGGGNGGCWGCAG
rev_primer=GGATTAGATACCCBDGTAGTC
front_primer_position=341f
rev_primer_position=805r
ref_db_trunc_len=465
metadata_file="filenames-biagi-raw-supercomp.tsv"
conda activate qiime2
# Import data
qiime tools import \
  --type  'SampleData[SequencesWithQuality]' \
  --input-path ${metadata_file} \
  --output-path ${author_name}-single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data ${author_name}-single-end-demux.qza \
  --o-visualization ${author_name}-single-end-demux.qzv
# Remove primers
qiime cutadapt trim-single \
  --i-demultiplexed-sequences ${author_name}-single-end-demux.qza \
  --p-cores 8 \
  --p-front ${front_primer} \
  --p-adapter ${rev_primer} \
  --o-trimmed-sequences ${author_name}-trimmed-single-end.qza
qiime demux summarize \
  --i-data ${author_name}-trimmed-single-end.qza \
  --o-visualization ${author_name}-trimmed-single-end.qzv
# denoising with dada2
# no truncation
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ${author_name}-trimmed-single-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences ${author_name}-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-table ${author_name}-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-denoising-stats ${author_name}-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza 
# visualize your dada2 stats
qiime metadata tabulate \
  --m-input-file ${author_name}-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization ${author_name}-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
# generate feature table (frequencies)
qiime feature-table summarize \
  --i-table ${author_name}-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization ${author_name}-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
# generate feature table (sequences)
qiime feature-table tabulate-seqs \
  --i-data ${author_name}-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization ${author_name}-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv 
# taxonomic classification with RESCRIPT
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138-1-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138-1-ssu-nr99-tax.qza
# RNA to DNA
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138-1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138-1-ssu-nr99-seqs.qza
qiime rescript cull-seqs \
  --i-sequences silva-138-1-ssu-nr99-seqs.qza \
  --o-clean-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza
# Extract reads
qiime feature-classifier extract-reads \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer ${front_primer} \
  --p-r-primer ${rev_primer} \
  --p-trunc-len ${ref_db_trunc_len} \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}.qza
# dereplicate extracted region  
qiime rescript dereplicate \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}.qza \
  --i-taxa silva-138-1-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}-uniq.qza \
  --o-dereplicated-taxa  silva-138-1-ssu-nr99-tax-${front_primer_position}-${rev_primer_position}-derep-uniq.qza
# train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${rev_primer_position}-uniq.qza \
  --i-reference-taxonomy silva-138-1-ssu-nr99-tax-${front_primer_position}-${rev_primer_position}-derep-uniq.qza \
  --o-classifier silva-138-1-ssu-nr99-classifier-${front_primer_position}-${rev_primer_position}-derep-uniq.qza
qiime feature-table filter-seqs \
  --i-data ${author_name}-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-metadata-file ${author_name}-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-where 'length(sequence) > 470' \
  --o-filtered-data ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime feature-table filter-features \
 --i-table ${author_name}-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
 --m-metadata-file ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
 --p-no-exclude-ids \
 --o-filtered-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime feature-table summarize \
  --i-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-table tabulate-seqs \
  --i-data ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
# test the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-1-ssu-nr99-classifier-${front_primer_position}-${rev_primer_position}-derep-uniq.qza \
  --i-reads ${author_name}-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-classification ${author_name}-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime taxa filter-table \
  --i-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --i-taxonomy ${author_name}-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza
qiime feature-table summarize \
  --i-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qzv
mv ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime feature-table summarize \
  --i-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-finalfiltercheck.qzv
qiime taxa filter-seqs \
  --i-sequences ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --i-taxonomy ${author_name}-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-sequences ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza
qiime feature-table tabulate-seqs \
  --i-data ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza \
  --o-visualization ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qzv
mv ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-nomitchlorarch.qza ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza
qiime feature-table tabulate-seqs \
  --i-data ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}-finalfiltercheck.qzv
qiime metadata tabulate \
  --m-input-file ${author_name}-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization ${author_name}-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime taxa barplot \
  --i-table ${author_name}-filtered-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --i-taxonomy ${author_name}-taxonomy-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --m-metadata-file ${metadata_file} \
  --o-visualization ${author_name}-taxa-bar-plots-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${author_name}-filtered-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-alignment ${author_name}-aligned-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-masked-alignment ${author_name}-masked-aligned-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-tree ${author_name}-unrooted-tree-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-rooted-tree ${author_name}-rooted-tree-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza