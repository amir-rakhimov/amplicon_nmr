#!/usr/bin/env bash
# 0. Install QIIME2
# wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
# mamba env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-linux-conda.yml
# 1. Set up the variables
time_var=$(date +%T |sed 's/:/_/g' )
date_var=$(date -I|sed 's/-//g')
date_time=${date_var}_${time_var}
### Truncation length for DADA2
trunc_len=234
project_home_dir=~/projects/amplicon_nmr
output_dir=${project_home_dir}/output/qiime/pooled-qiime
metadata_dir=${project_home_dir}/data/metadata/pooled-metadata
mkdir ${output_dir}/${date_time}-single-${trunc_len}
cd ${output_dir}/${date_time}-single-${trunc_len}
### Front and reverse primers
front_primer=CCTACGGGNGGCWGCAG
rev_primer=GACTACHVGGGTATCTAATCC
### Position of the front primer
front_primer_position=341f
### metadata file contains the fastq file paths
metadata_file=${metadata_dir}/filenames-single-pooled-raw-supercomp.tsv
### Truncation length for the taxonomic reference database
ref_db_trunc_len=234
conda activate  qiime2-amplicon-2024.2
# 2. Import the fastq files using the metadata file.
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
### The file shows minimum, maximum, median, mean, and total reads. 
### It also shows quality plots, frequency histograms, and number of reads in each sample
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
### Create a visualization to see that primers were removed: the sequences will get shorter.
### The contents are similar to the previous visualization
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
### Visualise the denoising-stats file: it shows input reads in each sample,
### how many reads passed the filtering, denoising, and chimera removal.
qiime metadata tabulate \
  --m-input-file pooled-single-stats-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-stats-trimmed-dada2-${trunc_len}.qzv
### Visualise the feature table (pooled-single-table): It shows the total number of samples, 
### number of features, and total frequency. You can check frequency per sample, too.
### With metadata, you can see the library size per metadata category (e.g. host, age, sex).
qiime feature-table summarize \
  --i-table pooled-single-table-trimmed-dada2-${trunc_len}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-table-trimmed-dada2-${trunc_len}.qzv
### Visualise representative-sequences (pooled-single-rep-seqs). 
### This will show you sequence lengths (min, max, mean, range, SD)
qiime feature-table tabulate-seqs \
  --i-data pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qzv
### 5. Create the reference database with RESCRIPT for taxonomic classification.
### Use SILVA v138.1
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
### 6. Fit a Naive Bayes classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-single-138-1-ssu-nr99-seqs-cleaned-${front_primer_position}-${ref_db_trunc_len}-uniq.qza \
  --i-reference-taxonomy silva-single-138-1-ssu-nr99-tax-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza \
  --o-classifier silva-single-138-1-ssu-nr99-classifier-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza
### 7.1 Filter representative sequences by length (remove shorter than 200).
### Metadata-based filtering: use the pooled-single-rep-seqs file as metadata because the list of 
### IDs to keep is determined based on metadata search  criteria rather than being provided by the user directly.
### This is achieved using the --p-where parameter in combination with the --m-metadata-file parameter. 
### The user provides a description of the samples that should be retained based on their metadata using 
### --p-where, where the syntax for this description is the SQLite WHERE-clause syntax.
### Result: pooled-single-filtered-rep-seqs
qiime feature-table filter-seqs \
  --i-data pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --m-metadata-file pooled-single-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --p-where 'length(sequence) > 199' \
  --o-filtered-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza
### 7.2 Filter the feature table (pooled-single-table becomes pooled-single-filtered-table).
### Metadata-based filtering: use FILTERED representative sequences (pooled-single-filtered-rep-seqs)
### as a list of metadata IDs. The --p-no-exclude-ids means features selected by metadata will be retained
qiime feature-table filter-features \
 --i-table pooled-single-table-trimmed-dada2-${trunc_len}.qza \
 --m-metadata-file pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
 --p-no-exclude-ids \
 --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza
### Visualise the filtered table (pooled-single-filtered-table).
### We see less features, lower frequency, etc
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qzv
### Visualise filtered representative sequences (pooled-single-filtered-rep-seqs)
### We see that sequence count decreased
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qzv
### 8. Classify filtered representative sequences (pooled-single-filtered-rep-seqs) with your classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva-single-138-1-ssu-nr99-classifier-${front_primer_position}-${ref_db_trunc_len}-derep-uniq.qza\
  --i-reads pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-classification pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza
### 9. Filter the filtered feature table (pooled-single-filtered-table) by taxonomy (pooled-single-taxonomy): 
### remove mitochondria, archaea, and chloroplast DNA. Use --i-taxonomy as filtering metadata.
### Output: pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch
qiime taxa filter-table \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --i-taxonomy pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qza
### Visualise the filtered table (pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch).
### We see less features, lower frequency, etc
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qzv
### Rename the new filtered table (pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch becomes
### pooled-single-filtered-table)
mv pooled-single-filtered-table-trimmed-dada2-${trunc_len}-nomitchlorarch.qza pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza
### Final visualization of the filtered feature table (pooled-single-filtered-table)
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-finalfiltercheck.qzv
### 10. Filter representative sequences (pooled-single-filtered-rep-seqs) by taxonomy (pooled-single-taxonomy):
### remove mitochondria, archaea, and chloroplast DNA. Use --i-taxonomy as filtering metadata.
### Output: pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch
qiime taxa filter-seqs \
  --i-sequences pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --i-taxonomy pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-sequences pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qza
### Visualise filtered representative sequences (pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch)
### We see that sequence count decreased
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qzv
### Rename the new filtered representative sequences (pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch becomes
### pooled-single-filtered-rep-seqs)
mv pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-nomitchlorarch.qza pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza
### Final visualization of the filtered representative sequences (pooled-single-filtered-rep-seqs)
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-finalfiltercheck.qzv
### Visualise taxonomy (not informative)
qiime metadata tabulate \
  --m-input-file pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --o-visualization pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qzv
### 11. Build a taxa barplot: use metadata file for visualization
qiime taxa barplot \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --i-taxonomy pooled-single-taxonomy-trimmed-dada2-${trunc_len}.qza \
  --m-metadata-file ${metadata_file} \
  --o-visualization pooled-single-taxa-bar-plots-trimmed-dada2-${trunc_len}.qzv
### 12. Align filtered representative sequences (pooled-single-filtered-rep-seqs) to a tree.
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-alignment pooled-single-aligned-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-masked-alignment pooled-single-masked-aligned-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --o-tree pooled-single-unrooted-tree-trimmed-dada2-${trunc_len}.qza \
  --o-rooted-tree pooled-single-rooted-tree-trimmed-dada2-${trunc_len}.qza