#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.2
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "$(date +"%F %H:%M:%S")"
for author in "${authors[@]}"; do
    echo "$(date +"%F %H:%M:%S")"
    echo "${author}"
    qiime2_output_dir="${OUTDIR}"/01-qiime2_output/"${author}"
    metadata_file="${metadata_dir}"/filenames-paired-"${author}".tsv
    mkdir -p "${qiime2_output_dir}"
    echo "$(date +"%F %H:%M:%S")"
    
    # 4. Denoise, dereplicate, and filter reads with DADA2.
    ### --denoise-single means we use single-end sequences
    ### --p-trim-left trims 5prime end 
    ### --p-trunc-len position at which sequences should be truncated. 
    ### It will truncate the 3prime end. Reads that are shorter than
    ### p-trunc-len will be discarded.
    echo "$(date +"%F %H:%M:%S")"
    # pairs=("270 210" "259 199" "284 203" "273 214" ) 
    if [[ "${author}" == "sibai" ]]; then
        pairs=("230 230" "230 225" "225 230" "225 225")
    else
        # pairs=("278 229" "278 225" "278 220" "278 215" "275 229" "275 225" "270 229" "270 225" "265 220" "230 230" "230 225" "225 230" "225 225" ) 
        pairs=("270 210" "259 199" "284 203" "273 214" ) 
    fi

    for pair in "${pairs[@]}"; do
        read -r test_trunc_len_f test_trunc_len_r <<< "${pair}"
        if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-table-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza ]]; then
            echo "Denoising, dereplicating, and filtering reads with DADA2 using coordinates: ${test_trunc_len_f} ${test_trunc_len_r}."
            echo "$(date +"%F %H:%M:%S")"
            qiime dada2 denoise-paired \
                --i-demultiplexed-seqs "${qiime2_output_dir}"/trimmed.qza \
                --p-n-threads "${nthreads_dada2}" \
                --verbose \
                --p-trim-left-f 0 \
                --p-trim-left-r 0 \
                --p-trunc-len-f "${test_trunc_len_f}" \
                --p-trunc-len-r "${test_trunc_len_r}" \
                --o-table "${qiime2_output_dir}"/trimmed-dada2-table-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza \
                --o-representative-sequences "${qiime2_output_dir}"/trimmed-dada2-rep-seqs-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza \
                --o-denoising-stats "${qiime2_output_dir}"/trimmed-dada2-stats-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza
        fi
        ### Visualise the denoising-stats file: it shows input reads in each sample,
        ### how many reads passed the filtering, denoising, and chimera removal.
        if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-stats-"${test_trunc_len_f}"_"${test_trunc_len_r}".qzv ]]; then
            echo "Visualising the denoised stats file using coordinates: ${test_trunc_len_f} ${test_trunc_len_r}"
            echo "$(date +"%F %H:%M:%S")"
            qiime metadata tabulate \
                --m-input-file "${qiime2_output_dir}"/trimmed-dada2-stats-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza \
                --o-visualization "${qiime2_output_dir}"/trimmed-dada2-stats-"${test_trunc_len_f}"_"${test_trunc_len_r}".qzv
        fi
        ### Visualise the feature table (table): It shows the total number of samples, 
        ### number of features, and total frequency. You can check frequency per sample, too.
        ### With metadata, you can see the library size per metadata category (e.g. host, age, sex).
        if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-table-"${test_trunc_len_f}"_"${test_trunc_len_r}".qzv ]]; then
            echo "$(date +"%F %H:%M:%S")"
            echo "Visualising the denoised feature table using coordinates: ${test_trunc_len_f} ${test_trunc_len_r}"
            qiime feature-table summarize \
                --i-table "${qiime2_output_dir}"/trimmed-dada2-table-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza \
                --m-sample-metadata-file "${metadata_file}" \
                --o-visualization "${qiime2_output_dir}"/trimmed-dada2-table-"${test_trunc_len_f}"_"${test_trunc_len_r}".qzv
        fi
        ### Visualise representative-sequences (rep-seqs). 
        ### This will show you sequence lengths (min, max, mean, range, SD)

        echo "$(date +"%F %H:%M:%S")"
        if [[ ! -f "${qiime2_output_dir}"/trimmed-dada2-rep-seqs-"${test_trunc_len_f}"_"${test_trunc_len_r}".qzv ]]; then
            echo "Visualising the denoised representative sequences"
            qiime feature-table tabulate-seqs \
                --i-data "${qiime2_output_dir}"/trimmed-dada2-rep-seqs-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza \
                --o-visualization "${qiime2_output_dir}"/trimmed-dada2-rep-seqs-"${test_trunc_len_f}"_"${test_trunc_len_r}".qzv
        fi
        

    done
    # Export statistics and representative sequences as tsv and fasta files
    ## Statistics
    mkdir -p "${qiime2_output_dir}"/exported-stats
    parallel  'unzip -j {1} "*/stats.tsv" -d "{2}/exported-stats/{1/.}" ' ::: "${qiime2_output_dir}"/*stats-*.qza ::: "${qiime2_output_dir}"
    ## Representative sequences
    mkdir -p "${qiime2_output_dir}"/exported-rep-seqs
    parallel  'unzip -j {1} "*/dna-sequences.fasta" -d "{2}/exported-rep-seqs/{1/.}" ' ::: "${qiime2_output_dir}"/*rep-seqs-*.qza ::: "${qiime2_output_dir}"
done

echo "$(date +"%F %H:%M:%S")"
