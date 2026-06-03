#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.2
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

echo "$(date +"%F %H:%M:%S")"
for author in "${authors[@]}"; do
    echo "$(date +"%F %H:%M:%S")"
    echo "${author}"
    qiime2_output_dir="${OUTDIR}"/01-qiime2_output/"${author}"
    metadata_file="${metadata_dir}"/filenames-paired-"${author}".tsv
    fastqc_raw_output_dir=data/1-qc/fastqc_raw/"${RUN_ID}"/"${author}"
    multiqc_raw_output_dir=data/1-qc/multiqc_raw/"${RUN_ID}"/"${author}"
    fastqc_trimmed_output_dir=data/1-qc/fastqc_trimmed/"${RUN_ID}"/"${author}"
    multiqc_trimmed_output_dir=data/1-qc/multiqc_trimmed/"${RUN_ID}"/"${author}"

    if [[ "${author}" == "yasuda" ]]; then
        fastq_dir=data/0-raw/fastq/for-ncbi/renamed
    else
        fastq_dir=data/0-raw/fastq/"${author}"-fastq/raw
    fi

    mkdir -p "${qiime2_output_dir}"
    mkdir -p "${fastqc_raw_output_dir}"
    mkdir -p "${multiqc_raw_output_dir}"
    mkdir -p "${fastqc_trimmed_output_dir}"
    mkdir -p "${multiqc_trimmed_output_dir}"
    echo "$(date +"%F %H:%M:%S")"

    # QC: Check sequence length distribution
    if [[ ! -f "${multiqc_raw_output_dir}"/R1_report.html ]]; then
        conda activate qc-tools 
        # Run FastQC and trim overrepresented sequences
        echo "Running FastQC on raw data"
        echo "$(date +"%F %H:%M:%S")"
        find "${fastq_dir}" -name '*.fastq.gz' | parallel -j "${nthreads_qc}" --verbose \
            fastqc {} --outdir "${fastqc_raw_output_dir}"  
        ## Run MultiQC
        echo "Running MultiQC on raw data"
        echo "$(date +"%F %H:%M:%S")"
        multiqc "${fastqc_raw_output_dir}"/ --outdir "${multiqc_raw_output_dir}" --filename combined_report
        multiqc "${fastqc_raw_output_dir}"/*R1_fastqc.zip --outdir "${multiqc_raw_output_dir}" --filename R1_report
        multiqc "${fastqc_raw_output_dir}"/*R2_fastqc.zip --outdir "${multiqc_raw_output_dir}" --filename R2_report
        conda deactivate
    fi

    # 2. Import the fastq files using the metadata file.
    if [[ ! -f "${qiime2_output_dir}"/demux.qza ]]; then
        echo "Importing FASTQ files"
        if [[ "${read_type}" == "single_end" ]]; then
            qiime tools import \
                --type  'SampleData[SequencesWithQuality]' \
                --input-path "${metadata_file}" \
                --output-path "${qiime2_output_dir}"/demux.qza \
                --input-format SingleEndFastqManifestPhred33V2 
        elif [[ "${read_type}" == "paired_end" ]]; then
            qiime tools import \
                --type  'SampleData[PairedEndSequencesWithQuality]' \
                --input-path "${metadata_file}" \
                --output-path "${qiime2_output_dir}"/demux.qza \
                --input-format PairedEndFastqManifestPhred33V2 
        fi 
    fi

    ### The command below will create a visualization file for quality check.
    ### The file shows minimum, maximum, median, mean, and total reads. 
    ### It also shows quality plots, frequency histograms, and number of reads in each sample
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${qiime2_output_dir}"/demux.qzv ]]; then
        echo "Visualising the raw data"
        qiime demux summarize \
            --i-data "${qiime2_output_dir}"/demux.qza \
            --o-visualization "${qiime2_output_dir}"/demux.qzv
    fi

    # 3. Remove the sequencing adapters with cutadapt.
    ### Specify the primer sequences. If you're using only front reads, use --p-front only
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${qiime2_output_dir}"/trimmed.qza ]]; then
        echo "Removing adapters with cutadapt"
        if [[ "${read_type}" == "single_end" ]]; then
            qiime cutadapt trim-single \
                --i-demultiplexed-sequences "${qiime2_output_dir}"/demux.qza \
                --p-discard-untrimmed True \
                --p-cores "${nthreads_qc}" \
                --p-front "${fwd_primer}" \
                --o-trimmed-sequences "${qiime2_output_dir}"/trimmed.qza
        elif [[ "${read_type}" == "paired_end" ]]; then
            qiime cutadapt trim-paired \
                --i-demultiplexed-sequences "${qiime2_output_dir}"/demux.qza \
                --p-cores "${nthreads_qc}" \
                --p-discard-untrimmed True \
                --p-front-f "${fwd_primer}" \
                --p-front-r "${rev_primer}" \
                --o-trimmed-sequences "${qiime2_output_dir}"/trimmed.qza
        fi
    fi

    ### Create a visualization to see that primers were removed: the sequences will get shorter.
    ### The contents are similar to the previous visualization
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${qiime2_output_dir}"/trimmed.qzv ]]; then
        echo "Visualising the trimmed data"
        qiime demux summarize \
            --i-data "${qiime2_output_dir}"/trimmed.qza \
            --o-visualization "${qiime2_output_dir}"/trimmed.qzv
    fi

    
    if [[ ! -d "${qiime2_output_dir}"/exported-trimmed-reads ]]; then
        qiime tools export \
            --input-path "${qiime2_output_dir}"/trimmed.qza \
            --output-path "${qiime2_output_dir}"/exported-trimmed-reads
    fi
    # Run FastQC to check the trimming result

    if [[ ! -f "${multiqc_trimmed_output_dir}"/R1_report.html ]]; then
        conda activate qc-tools 
        # Run FastQC to check the trimming result
        echo "Running FastQC on trimmed data"
        echo "$(date +"%F %H:%M:%S")"
        find "${qiime2_output_dir}"/exported-trimmed-reads -name '*.fastq.gz' | parallel -j "${nthreads_qc}" --verbose \
            fastqc {} --outdir "${fastqc_trimmed_output_dir}"  
        ## Run MultiQC
        echo "Running MultiQC on trimmed data"
        echo "$(date +"%F %H:%M:%S")"
        # multiqc "${fastqc_trimmed_output_dir}"/ --outdir "${multiqc_trimmed_output_dir}" --filename combined_report
        multiqc "${fastqc_trimmed_output_dir}"/*_R1_001_fastqc.zip --outdir "${multiqc_trimmed_output_dir}" --filename R1_report
        multiqc "${fastqc_trimmed_output_dir}"/*_R2_001_fastqc.zip --outdir "${multiqc_trimmed_output_dir}" --filename R2_report
        conda deactivate
    fi

done
echo "$(date +"%F %H:%M:%S")"
