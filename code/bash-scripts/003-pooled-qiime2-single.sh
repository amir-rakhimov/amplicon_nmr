#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 16
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20260209_16-32-pooled-qiime2-single
#SBATCH --output jobreports/20260209_16-32-pooled-qiime2-single-out.txt
#SBATCH --error jobreports/20260209_16-32-pooled-qiime2-single-out.txt
# change the shebang to /usr/bin/env bash if needed
# Use supercomputer because fit-classifier-naive-bayes and classify-sklearn were killed on my desktop
# 0. Install QIIME2
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# wget https://data.qiime2.org/distro/core/qiime2-2022.8-py38-linux-conda.yml -O ~/qiime2-2022.8-py38-linux-conda.yml
# sed "s/^  //" ~/qiime2-2022.8-py38-linux-conda.yml > ~/qiime2-2022.8-py38-linux-conda-edited.yml
# mamba create -yn qiime2-2022.8 --file ~/qiime2-2022.8-py38-linux-conda-edited.yml

# This one worked:
# wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml -O ~/qiime2-amplicon-2024.2-py38-linux-conda.yml
# mamba env create -n qiime2-amplicon-2024.2 --file ~/qiime2-amplicon-2024.2-py38-linux-conda.yml
# 1. Set up the variables
time_var=$(date +%T |sed 's/:/_/g' )
date_var=$(date -I|sed 's/-//g')
date_time=${date_var}_${time_var}
# date_time=20240425_02_57_13
# Start of the script execution
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
### Truncation length for DADA2
trunc_len=234
project_home_dir=~/projects/amplicon_nmr
output_dir="${project_home_dir}"/output/qiime/pooled-qiime
metadata_dir="${project_home_dir}"/data/metadata/pooled-metadata
nthreads=14
mkdir -p "${output_dir}"/"${date_time}"-single-"${trunc_len}"
qiime_output_dir="${output_dir}"/"${date_time}"-single-"${trunc_len}"
### Front and reverse primers
fwd_primer=CCTACGGGNGGCWGCAG
rev_primer=GACTACHVGGGTATCTAATCC
### Position of the front primer
fwd_primer_position=341f
### metadata file contains the fastq file paths
metadata_file="${metadata_dir}"/filenames-single-pooled-raw-supercomp.tsv
### Truncation length for the taxonomic reference database
ref_db_trunc_len=234
conda activate qiime2-amplicon-2024.2
# conda activate qiime2-2022.8
# 2. Import the fastq files using the metadata file.
### We need to specify the input format. Since we're using only
### front reads (single reads), we need to specify SingleEndFastqManifestPhred33V2
echo "Importing FASTQ files"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime tools import \
  --type  'SampleData[SequencesWithQuality]' \
  --input-path "${metadata_file}" \
  --output-path "${qiime_output_dir}"/"${date_time}"-pooled-single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 
### The command below will create a visualization file for quality check.
### Normally, you should open the file and decide the truncation length. 
### But I have already explored the final output with different truncation lentgths,
### so it's already set at the beginning of this script. 
### The file shows minimum, maximum, median, mean, and total reads. 
### It also shows quality plots, frequency histograms, and number of reads in each sample
echo "Visualising the raw data"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime demux summarize \
  --i-data "${qiime_output_dir}"/"${date_time}"-pooled-single-end-demux.qza \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-end-demux.qzv
# 3. Remove the sequencing adapters with cutadapt.
### Specify the primer sequences. We're using only front reads, so we specify
### only --p-front
echo "Removing adapters with cutadapt"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime cutadapt trim-single \
  --i-demultiplexed-sequences "${qiime_output_dir}"/"${date_time}"-pooled-single-end-demux.qza \
  --p-cores "${nthreads}" \
  --p-front "${fwd_primer}" \
  --o-trimmed-sequences "${qiime_output_dir}"/"${date_time}"-pooled-single-end-trimmed.qza
### Create a visualization to see that primers were removed: the sequences will get shorter.
### The contents are similar to the previous visualization
echo "Visualising the trimmed data"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime demux summarize \
  --i-data "${qiime_output_dir}"/"${date_time}"-pooled-single-end-trimmed.qza \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-end-trimmed.qzv
# 4. Denoise, dereplicate, and filter reads with DADA2.
### --denoise-single means we use single-end sequences
### --p-trim-left trims 5prime end 
### --p-trunc-len position at which sequences should be truncated. 
### It will truncate the 3prime end. Reads that are shorter than
### p-trunc-len will be discarded.
echo "Denoising, dereplicating, and filtering reads with DADA2."
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime dada2 denoise-single \
  --p-n-threads "${nthreads}" \
  --i-demultiplexed-seqs "${qiime_output_dir}"/"${date_time}"-pooled-single-end-trimmed.qza \
  --p-trim-left 0 \
  --p-trunc-len "${trunc_len}" \
  --o-table "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}".qza \
  --o-representative-sequences "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}".qza \
  --o-denoising-stats "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-stats-"${trunc_len}".qza
### Visualise the denoising-stats file: it shows input reads in each sample,
### how many reads passed the filtering, denoising, and chimera removal.
echo "Visualising the denoised stats file"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime metadata tabulate \
  --m-input-file "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-stats-"${trunc_len}".qza \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-stats-"${trunc_len}".qzv
### Visualise the feature table (pooled-single-table): It shows the total number of samples, 
### number of features, and total frequency. You can check frequency per sample, too.
### With metadata, you can see the library size per metadata category (e.g. host, age, sex).
echo "Visualising the denoised feature table"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-table summarize \
  --i-table "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}".qza \
  --m-sample-metadata-file "${metadata_file}" \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}".qzv
### Visualise representative-sequences (pooled-single-rep-seqs). 
### This will show you sequence lengths (min, max, mean, range, SD)
echo "Visualising the denoised representative sequences"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-table tabulate-seqs \
  --i-data "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}".qza \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}".qzv
### 5. Create the reference database with RESCRIPT for taxonomic classification.
### Use SILVA v138.1
echo "Creating the reference database with RESCRIPT for taxonomic classification"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime rescript get-silva-data \
  --p-version '138.1' \
  --p-target 'SSURef_NR99' \
  --p-include-species-labels \
  --o-silva-sequences "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-rna-seqs.qza \
  --o-silva-taxonomy "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-tax.qza
echo "Reverse transcribing RNA to DNA"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime rescript reverse-transcribe \
  --i-rna-sequences "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-seqs.qza
echo "Cleaning sequences: remove sequences that contain 5 or more ambiguous bases \
  (IUPAC compliant ambiguity bases) \
  and any homopolymers that are 8 or more bases in length. \
  These are the default parameters. See the --help text for more details."
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime rescript cull-seqs \
  --i-sequences "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-seqs.qza \
  --o-clean-sequences "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-seqs-cleaned.qza
echo "# Extracting reads"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-classifier extract-reads \
  --i-sequences "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer "${fwd_primer}" \
  --p-r-primer "${rev_primer}" \
  --p-trunc-len "${ref_db_trunc_len}" \
  --p-n-jobs "${nthreads}" \
  --p-read-orientation 'forward' \
  --o-reads "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-seqs-cleaned-"${fwd_primer_position}"-"${ref_db_trunc_len}".qza
echo "Dereplicating the extracted region"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime rescript dereplicate \
  --i-sequences "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-seqs-cleaned-"${fwd_primer_position}"-"${ref_db_trunc_len}".qza \
  --i-taxa "${qiime_output_dir}"/"${date_time}"-silva-138-1-ssu-nr99-tax.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-seqs-cleaned-"${fwd_primer_position}"-"${ref_db_trunc_len}"-uniq.qza \
  --o-dereplicated-taxa "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-tax-"${fwd_primer_position}"-"${ref_db_trunc_len}"-derep-uniq.qza
### 6. Fit a Naive Bayes classifier
echo "Training a Naive Bayes classifier"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-seqs-cleaned-"${fwd_primer_position}"-"${ref_db_trunc_len}"-uniq.qza \
  --i-reference-taxonomy "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-tax-"${fwd_primer_position}"-"${ref_db_trunc_len}"-derep-uniq.qza \
  --o-classifier "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-classifier-"${fwd_primer_position}"-"${ref_db_trunc_len}"-derep-uniq.qza
### 7.1 Filter representative sequences by length (remove shorter than 200).
### Metadata-based filtering: use the pooled-single-rep-seqs file as metadata because the list of 
### IDs to keep is determined based on metadata search  criteria rather than being provided by the user directly.
### This is achieved using the --p-where parameter in combination with the --m-metadata-file parameter. 
### The user provides a description of the samples that should be retained based on their metadata using 
### --p-where, where the syntax for this description is the SQLite WHERE-clause syntax.
### Result: pooled-single-rep-seqs-filtered
echo "Filtering the representative sequences by length (remove shorter than 200 bp)"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-table filter-seqs \
  --i-data "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}".qza \
  --m-metadata-file "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}".qza \
  --p-where 'length(sequence) > 199' \
  --o-filtered-data "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered.qza
### 7.2 Filter the feature table (pooled-single-table becomes pooled-single-table-filtered).
### Metadata-based filtering: use FILTERED representative sequences (pooled-single-rep-seqs-filtered)
### as a list of metadata IDs. The --p-no-exclude-ids means features selected by metadata will be retained
echo "Filtering the feature table using FILTERED representative sequences as metadata"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-table filter-features \
 --i-table "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}".qza \
 --m-metadata-file "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered.qza \
 --p-no-exclude-ids \
 --o-filtered-table "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}"-filtered.qza
echo "Visualising the filtered feature table"
### Visualise the filtered table (pooled-single-filtered-table).
### We see less features, lower frequency, etc
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-table summarize \
  --i-table "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}"-filtered.qza \
  --m-sample-metadata-file "${metadata_file}" \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}"-filtered.qzv
echo "Visualising the filtered representative sequences"
### Visualise filtered representative sequences (pooled-single-filtered-rep-seqs)
### We see that sequence count decreased
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-table tabulate-seqs \
  --i-data "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered.qza \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered.qzv
### 8. Classify filtered representative sequences (pooled-single-filtered-rep-seqs) with your classifier
echo "Classifying the filtered representative sequences"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime feature-classifier classify-sklearn \
  --i-classifier "${qiime_output_dir}"/"${date_time}"-silva-single-138-1-ssu-nr99-classifier-"${fwd_primer_position}"-"${ref_db_trunc_len}"-derep-uniq.qza \
  --i-reads "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered.qza \
  --o-classification "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-"${trunc_len}"-filtered-taxonomy.qza
### 9. Build a taxa barplot: use metadata file for visualization
echo "Visualising the taxonomy with a barplot"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime taxa barplot \
  --i-table "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-table-"${trunc_len}"-filtered.qza \
  --i-taxonomy "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-"${trunc_len}"-filtered-taxonomy.qza \
  --m-metadata-file "${metadata_file}" \
  --o-visualization "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-"${trunc_len}"-taxa-bar-plots.qzv
### 10. Align filtered representative sequences (pooled-single-filtered-rep-seqs) to a tree.
echo "Aligning the filtered representative sequences to a tree"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered.qza \
  --o-alignment "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered-aligned.qza \
  --o-masked-alignment "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-rep-seqs-"${trunc_len}"-filtered-aligned-masked.qza \
  --o-tree "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-"${trunc_len}"-filtered-unrooted-tree.qza \
  --o-rooted-tree "${qiime_output_dir}"/"${date_time}"-pooled-single-trimmed-dada2-"${trunc_len}"-filtered-rooted-tree.qza
# End of the script execution
echo "End"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
