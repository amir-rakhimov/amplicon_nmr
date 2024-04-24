#!/usr/bin/env bash
# 1. Set up the variables
### Truncation length for DADA2
trunc_len=234
### Front and reverse primers
front_primer=CCTACGGGNGGCWGCAG
rev_primer=GACTACHVGGGTATCTAATCC
### Position of the front primer
front_primer_position=341f
### metadata file contains the fastq file paths
metadata_file="filenames-single-pooled-raw-supercomp.tsv"
### Truncation length for the taxonomic reference database
ref_db_trunc_len=234
conda activate qiime2
# 2 Import the fastq files using the metadata file.
### We need to specify the input format. Since we're using only
### front reads (single reads), we need to specify SingleEndFastqManifestPhred33V2
qiime tools import \
  --type  'SampleData[SequencesWithQuality]' \
  --input-path ${metadata_file} \
  --output-path pooled-single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 
### The command below will create a visualization file for quality check.
### Normally, you should open the file and decide the truncation length. 
### But I have already explored the final output with different truncation lentgths,
### so it's already set at the beginning of this script. 
qiime demux summarize \
  --i-data pooled-single-end-demux.qza \
  --o-visualization pooled-single-end-demux.qzv
# 3. Remove the sequencing adapters with cutadapt.
### Specify the primer sequences. We're using only front reads, so we specify
### only --p-front
qiime cutadapt trim-single \
  --i-demultiplexed-sequences pooled-single-end-demux.qza \
  --p-cores 8 \
  --p-front ${front_primer} \
  --o-trimmed-sequences pooled-trimmed-single-end.qza
### Create a visualization to see that primers were removed: the sequences will get shorter
qiime demux summarize \
  --i-data pooled-trimmed-single-end.qza \
  --o-visualization pooled-trimmed-single-end.qzv
# 4. Denoise, dereplicate, and filter reads with DADA2.
### --denoise-single means we use single-end sequences
### --p-trim-left trims 5prime end 
### --p-trunc-len position at which sequences should be truncated. 
### It will truncate the 3prime end. Reads that are shorter than
### p-trunc-len will be discarded.
qiime dada2 denoise-single \
  --i-demultiplexed-seqs pooled-trimmed-single-end.qza \
  --p-trim-left 0 \
  --p-trunc-len ${trunc_len} \
  --o-table pooled-single-table-trimmed-dada2-${trunc_len}.qza \
  --o-representative-sequences pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-denoising-stats pooled-single-stats-trimmed-dada2-${trunc_len}.qza
### Visualise the denoising-stats file
qiime metadata tabulate \
  --m-input-file pooled-single-stats-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-stats-trimmed-dada2-${trunc_len}.qzv
### Visualise the feature table. Use the metadata file to check library size by metadata
qiime feature-table summarize \
  --i-table pooled-single-table-trimmed-dada2-${trunc_len}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-table-trimmed-dada2-${trunc_len}.qzv
### Visualise representative-sequences. This will show you sequence lengths
qiime feature-table tabulate-seqs \
  --i-data pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qzv
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
  --o-reads silva-single-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${ref_db_trunc_len}.qza
qiime rescript dereplicate \
  --i-sequences silva-single-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${ref_db_trunc_len}.qza \
  --i-taxa silva-138-1-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-single-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${ref_db_trunc_len}-uniq.qza \
  --o-dereplicated-taxa  silva-single-138-1-ssu-nr99-tax-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-single-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${ref_db_trunc_len}-uniq.qza \
  --i-reference-taxonomy silva-single-138-1-ssu-nr99-tax-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza \
  --o-classifier silva-single-138-1-ssu-nr99-classifier-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza
qiime feature-table filter-seqs \
  --i-data pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --m-metadata-file pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --p-where 'length(sequence) > 199' \
  --o-filtered-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza
qiime feature-table filter-features \
 --i-table pooled-single-table-trimmed-dada2-${trunc_len}.qza \
 --m-metadata-file pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
 --p-no-exclude-ids \
 --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qzv
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier silva-single-138-1-ssu-nr99-classifier-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza\
  --i-reads pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-classification pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza
qiime taxa filter-table \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --i-taxonomy pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qza
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qzv
mv pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qza pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-finalfiltercheck.qzv
qiime taxa filter-seqs \
  --i-sequences pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --i-taxonomy pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-sequences pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qza
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qzv
mv pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qza pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-finalfiltercheck.qzv
qiime metadata tabulate \
  --m-input-file pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qzv
qiime taxa barplot \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --i-taxonomy pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --m-metadata-file ${metadata_file} \
  --o-visualization pooled-single-taxa-bar-plots-trimmed-dada2-${trunc_len}.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-alignment pooled-single-aligned-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-masked-alignment pooled-single-masked-aligned-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-tree pooled-single-unrooted-tree-trimmed-dada2-${trunc_len}.qza \
  --o-rooted-tree pooled-single-rooted-tree-trimmed-dada2-${trunc_len}.qza