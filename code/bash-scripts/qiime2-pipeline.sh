#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
# 0. Install QIIME2
echo "$(date +"%F %H:%M:%S")"
# wget https://data.qiime2.org/distro/core/qiime2-2022.8-py38-linux-conda.yml -O $HOME/qiime2-2022.8-py38-linux-conda.yml
# sed "s/^  //" $HOME/qiime2-2022.8-py38-linux-conda.yml > $HOME/qiime2-2022.8-py38-linux-conda-edited.yml
# mamba create -yn qiime2-2022.8 --file $HOME/qiime2-2022.8-py38-linux-conda-edited.yml

# This one worked:
# wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml -O $HOME/qiime2-amplicon-2024.2-py38-linux-conda.yml
# mamba env create -n qiime2-amplicon-2024.2 --file $HOME/qiime2-amplicon-2024.2-py38-linux-conda.yml

# QC: Check sequence length distribution
# if [[ "${qc_only}" == "TRUE" ]] ; then
#     conda activate qc-tools 
#     # Run FastQC and trim overrepresented sequences
#     mkdir -p "${fastqc_raw_output_dir}"
#     mkdir -p "${multiqc_raw_output_dir}"
#     echo "Running FastQC on raw data"
#     intermediate_date_time=$(date +"%F %H:%M:%S")
#     echo "${intermediate_date_time}"
#     find "${fastq_dir}" -name '*.fastq.gz' | grep -v "sibai" | parallel -j "${nthreads}" --verbose \
#         fastqc {} --outdir "${fastqc_raw_output_dir}"  
#     ## Run MultiQC
#     echo "Running MultiQC on raw data"
#     intermediate_date_time=$(date +"%F %H:%M:%S")
#     echo "${intermediate_date_time}"
#     multiqc "${fastqc_raw_output_dir}"/ --outdir "${multiqc_raw_output_dir}" --filename combined_report
#     multiqc "${fastqc_raw_output_dir}"/*R1_fastqc.zip --outdir "${multiqc_raw_output_dir}" --filename R1_report
#     multiqc "${fastqc_raw_output_dir}"/*R2_fastqc.zip --outdir "${multiqc_raw_output_dir}" --filename R2_report
#     exit 1
# fi

# # 1. Run QIIME2
# echo "$(date +"%F %H:%M:%S")"
# conda activate qiime2-amplicon-2024.2
# # conda activate qiime2-2022.8
# # 2. Import the fastq files using the metadata file.
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/demux.qza ]]; then
#     echo "Importing FASTQ files"
#     if [[ "${read_type}" == "single_end" ]]; then
#         qiime tools import \
#             --type  'SampleData[SequencesWithQuality]' \
#             --input-path "${metadata_file}" \
#             --output-path "${qiime2_output_dir}"/demux.qza \
#             --input-format SingleEndFastqManifestPhred33V2 
#     elif [[ "${read_type}" == "paired_end" ]]; then
#         qiime tools import \
#             --type  'SampleData[PairedEndSequencesWithQuality]' \
#             --input-path "${metadata_file}" \
#             --output-path "${qiime2_output_dir}"/demux.qza \
#             --input-format PairedEndFastqManifestPhred33V2 
#     fi 
# fi

# ### The command below will create a visualization file for quality check.
# ### The file shows minimum, maximum, median, mean, and total reads. 
# ### It also shows quality plots, frequency histograms, and number of reads in each sample
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/demux.qzv ]]; then
#     echo "Visualising the raw data"
#     qiime demux summarize \
#         --i-data "${qiime2_output_dir}"/demux.qza \
#         --o-visualization "${qiime2_output_dir}"/demux.qzv
# fi

# # 3. Remove the sequencing adapters with cutadapt.
# ### Specify the primer sequences. If you're using only front reads, use --p-front only
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/trimmed.qza ]]; then
#     echo "Removing adapters with cutadapt"
#     if [[ "${read_type}" == "single_end" ]]; then
#         qiime cutadapt trim-single \
#             --i-demultiplexed-sequences "${qiime2_output_dir}"/demux.qza \
#             --p-cores "${nthreads}" \
#             --p-front "${fwd_primer}" \
#             --o-trimmed-sequences "${qiime2_output_dir}"/trimmed.qza
#     elif [[ "${read_type}" == "paired_end" ]]; then
#         qiime cutadapt trim-paired \
#             --i-demultiplexed-sequences "${qiime2_output_dir}"/demux.qza \
#             --p-cores "${nthreads}" \
#             --p-discard-untrimmed True \
#             --p-front-f "${fwd_primer}" \
#             --p-front-r "${rev_primer}" \
#             --o-trimmed-sequences "${qiime2_output_dir}"/trimmed.qza
#     fi
# fi

# ### Create a visualization to see that primers were removed: the sequences will get shorter.
# ### The contents are similar to the previous visualization
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/trimmed.qzv ]]; then
#     echo "Visualising the trimmed data"
#     qiime demux summarize \
#         --i-data "${qiime2_output_dir}"/trimmed.qza \
#         --o-visualization "${qiime2_output_dir}"/trimmed.qzv
# fi

# qiime tools export \
#     --input-path "${qiime2_output_dir}"/trimmed.qza \
#     --output-path "${qiime2_output_dir}"/exported-trimmed-reads
# conda activate qc-tools 
# # Run FastQC to check the trimming result
# mkdir -p "${fastqc_trimmed_output_dir}"
# mkdir -p "${multiqc_trimmed_output_dir}"
# echo "Running FastQC on trimmed data"
# echo "$(date +"%F %H:%M:%S")"
# find "${qiime2_output_dir}"/exported-trimmed-reads -name '*.fastq.gz' | parallel -j "${nthreads}" --verbose \
#     fastqc {} --outdir "${fastqc_trimmed_output_dir}"  
# ## Run MultiQC
# echo "Running MultiQC on trimmed data"
# echo "$(date +"%F %H:%M:%S")"
# # multiqc "${fastqc_trimmed_output_dir}"/ --outdir "${multiqc_trimmed_output_dir}" --filename combined_report
# multiqc "${fastqc_trimmed_output_dir}"/*_R1_001_fastqc.zip --outdir "${multiqc_trimmed_output_dir}" --filename R1_report
# multiqc "${fastqc_trimmed_output_dir}"/*_R2_001_fastqc.zip --outdir "${multiqc_trimmed_output_dir}" --filename R2_report
# conda deactivate

# 4. Denoise, dereplicate, and filter reads with DADA2.
### --denoise-single means we use single-end sequences
### --p-trim-left trims 5prime end 
### --p-trunc-len position at which sequences should be truncated. 
### It will truncate the 3prime end. Reads that are shorter than
### p-trunc-len will be discarded.
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-table.qza ]]; then
#     echo "Denoising, dereplicating, and filtering reads with DADA2."
#     if [[ "${read_type}" == "single_end" ]]; then
#         qiime dada2 denoise-single \
#             --p-n-threads "${nthreads_dada2}" \
#             --i-demultiplexed-seqs "${qiime2_output_dir}"/trimmed.qza \
#             --p-trim-left 0 \
#             --p-trunc-len "${trunc_len_f}" \
#             --o-table "${qiime2_output_dir}"/trimmed-dada2-table.qza \
#             --o-representative-sequences "${qiime2_output_dir}"/trimmed-dada2-rep-seqs.qza \
#             --o-denoising-stats "${qiime2_output_dir}"/trimmed-dada2-stats.qza
#     elif [[ "${read_type}" == "paired_end" ]]; then
#         qiime dada2 denoise-paired \
#             --i-demultiplexed-seqs "${qiime2_output_dir}"/trimmed.qza \
#             --p-n-threads "${nthreads_dada2}" \
#             --p-trim-left-f 0 \
#             --p-trim-left-r 0 \
#             --p-trunc-len-f "${trunc_len_f}" \
#             --p-trunc-len-r "${trunc_len_r}" \
#             --o-table "${qiime2_output_dir}"/trimmed-dada2-table.qza \
#             --o-representative-sequences "${qiime2_output_dir}"/trimmed-dada2-rep-seqs.qza \
#             --o-denoising-stats "${qiime2_output_dir}"/trimmed-dada2-stats.qza
#     fi
# fi

# ### Visualise the denoising-stats file: it shows input reads in each sample,
# ### how many reads passed the filtering, denoising, and chimera removal.
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-stats.qzv ]]; then
#     echo "Visualising the denoised stats file"
#     qiime metadata tabulate \
#         --m-input-file "${qiime2_output_dir}"/trimmed-dada2-stats.qza \
#         --o-visualization "${qiime2_output_dir}"/trimmed-dada2-stats.qzv
# fi

# ### Visualise the feature table (table): It shows the total number of samples, 
# ### number of features, and total frequency. You can check frequency per sample, too.
# ### With metadata, you can see the library size per metadata category (e.g. host, age, sex).
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-table.qzv ]]; then
#     echo "Visualising the denoised feature table"
#     qiime feature-table summarize \
#         --i-table "${qiime2_output_dir}"/trimmed-dada2-table.qza \
#         --m-sample-metadata-file "${metadata_file}" \
#         --o-visualization "${qiime2_output_dir}"/trimmed-dada2-table.qzv
# fi

# ### Visualise representative-sequences (rep-seqs). 
# ### This will show you sequence lengths (min, max, mean, range, SD)
# echo "$(date +"%F %H:%M:%S")"
# if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-rep-seqs.qzv ]]; then
#     echo "Visualising the denoised representative sequences"
#     qiime feature-table tabulate-seqs \
#         --i-data "${qiime2_output_dir}"/trimmed-dada2-rep-seqs.qza \
#         --o-visualization "${qiime2_output_dir}"/trimmed-dada2-rep-seqs.qzv
# fi
