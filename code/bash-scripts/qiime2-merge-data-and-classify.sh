#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.2
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "$(date +"%F %H:%M:%S")"
# pairs=("278 229" "278 225" "278 220" "278 215" "275 229" "275 225" "270 229" "270 225" "265 220" "230 230" "230 225" "225 230" "225 225" "270 210" "259 199" "284 203" "273 214" ) 
naive_bayes_classifier="${silva_data_dir}"/naive-bayes-classifier.qza
# final_pairs=("230 225" "225 225") 
mkdir -p "${qiime2_output_dir}"/merged-data
for pair in "${final_pairs[@]}"; do
    read -r test_trunc_len_f test_trunc_len_r <<< "${pair}"
    echo processing pair "${test_trunc_len_f}" "${test_trunc_len_r}"
    merged_data_dir="${qiime2_output_dir}"/merged-data/"${test_trunc_len_f}"_"${test_trunc_len_r}"
    mkdir -p "${merged_data_dir}"
    # Create an array of table paths
    mapfile -t final_tables < <(
        ls "${qiime2_output_dir}"/*/trimmed-dada2-table-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza |\
                        grep -E "$(export IFS="|"; echo "${final_authors[*]}")" 
    )
    # Build --i-tables args dynamically
    merge_tables_args=()
    for t in "${final_tables[@]}"; do
        merge_tables_args+=(--i-tables "$t")
    done

    # Create an array of rep-seqs paths
    mapfile -t final_rep_seqs < <(
        ls "${qiime2_output_dir}"/*/trimmed-dada2-rep-seqs-"${test_trunc_len_f}"_"${test_trunc_len_r}".qza |\
                        grep -E "$(export IFS="|"; echo "${final_authors[*]}")" 
    )
    # Build --i-data dynamically
    merge_rep_seqs_args=()
    for t in "${final_rep_seqs[@]}"; do
        merge_rep_seqs_args+=(--i-data "$t")
    done    
    # Define paths
    merged_table="${merged_data_dir}"/trimmed-dada2-table-merged.qza 
    merged_rep_seqs="${merged_data_dir}"/trimmed-dada2-rep-seqs-merged.qza

    merged_rep_seqs_visualisation=$(echo "${merged_rep_seqs}" | sed "s/qza/qzv/")

    merged_filtered_table="${merged_data_dir}"/trimmed-dada2-table-merged-filtered.qza
    merged_filtered_table_visualisation=$(echo "${merged_filtered_table}" | sed "s/qza/qzv/")
    
    merged_filtered_rep_seqs="${merged_data_dir}"/trimmed-dada2-rep-seqs-merged-filtered.qza
    merged_filtered_rep_seqs_visualisation=$(echo "${merged_filtered_rep_seqs}" | sed "s/qza/qzv/")

    taxonomy="${merged_data_dir}"/trimmed-dada2-merged-taxonomy.qza
    alignment="${merged_data_dir}"/trimmed-dada2-rep-seqs-merged-filtered-aligned.qza
    masked_alignment="${merged_data_dir}"/trimmed-dada2-rep-seqs-merged-filtered-aligned-masked.qza
    unrooted_tree="${merged_data_dir}"/trimmed-dada2-merged-filtered-unrooted-tree.qza
    rooted_tree="${merged_data_dir}"/trimmed-dada2-merged-filtered-rooted-tree.qza

    if [[ ! -f "${merged_table}" ]]; then
        echo "Merging tables with coordinates: ${test_trunc_len_f} ${test_trunc_len_r}."
        qiime feature-table merge \
            --i-tables "${merge_tables_args[@]}" \
            --o-merged-table "${merged_table}"
    fi

    if [[ ! -f "${merged_rep_seqs}" ]]; then
        echo "Merging representative sequences with coordinates: ${test_trunc_len_f} ${test_trunc_len_r}."
        qiime feature-table merge-seqs \
            --i-data "${merge_rep_seqs_args[@]}" \
            --o-merged-data "${merged_rep_seqs}"
    fi

    if [[ ! -f "${merged_rep_seqs_visualisation}" ]]; then
        echo "Visualising the merged representative sequences with coordinates: ${test_trunc_len_f} ${test_trunc_len_r}."
        qiime feature-table tabulate-seqs \
            --i-data "${merged_rep_seqs}" \
            --o-visualization "${merged_rep_seqs_visualisation}"
    fi

    ### 7.1 Filter representative sequences by length.
    ### Metadata-based filtering: use the rep-seqs file as metadata because the list of 
    ### IDs to keep is determined based on metadata search  criteria rather than being provided by the user directly.
    ### This is achieved using the --p-where parameter in combination with the --m-metadata-file parameter. 
    ### The user provides a description of the samples that should be retained based on their metadata using 
    ### --p-where, where the syntax for this description is the SQLite WHERE-clause syntax.
    ### Result: rep-seqs-filtered
    # echo "Filtering the representative sequences by length"
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${merged_filtered_rep_seqs}" ]]; then
        echo "Filtering the representative sequences by length (remove shorter than 380 and longer than 440 bp)"
        qiime feature-table filter-seqs \
            --i-data "${merged_rep_seqs}" \
            --m-metadata-file "${merged_rep_seqs}" \
            --p-where 'length(sequence) > 380 AND length(sequence)<440' \
            --o-filtered-data "${merged_filtered_rep_seqs}"
    fi
    ### 7.2 Filter the feature table (table becomes table-filtered).
    ### Metadata-based filtering: use FILTERED representative sequences (rep-seqs-filtered)
    ### as a list of metadata IDs. The --p-no-exclude-ids means features selected by metadata will be retained
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${merged_filtered_table}" ]]; then
        echo "Filtering the feature table using FILTERED representative sequences as metadata"
        qiime feature-table filter-features \
            --i-table "${merged_table}" \
            --m-metadata-file "${merged_filtered_rep_seqs}" \
            --p-no-exclude-ids \
            --o-filtered-table "${merged_filtered_table}"
    fi
    ### Visualise the filtered table (filtered-table).
    ### We see less features, lower frequency, etc
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${merged_filtered_table_visualisation}" ]]; then
        echo "Visualising the filtered feature table"
        qiime feature-table summarize \
            --i-table "${merged_filtered_table}" \
            --m-sample-metadata-file "${metadata_file}" \
            --o-visualization "${merged_filtered_table_visualisation}"
    fi
    ### Visualise filtered representative sequences (filtered-rep-seqs)
    ### We see that sequence count decreased
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${merged_filtered_rep_seqs_visualisation}" ]]; then
        echo "Visualising the filtered representative sequences"
        qiime feature-table tabulate-seqs \
            --i-data "${merged_filtered_rep_seqs}" \
            --o-visualization "${merged_filtered_rep_seqs_visualisation}"
    fi
    ### 8. Classify filtered representative sequences (filtered-rep-seqs) with your classifier
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${taxonomy}" ]]; then
        echo "Classifying the filtered representative sequences"
        qiime feature-classifier classify-sklearn \
            --i-classifier "${naive_bayes_classifier}" \
            --i-reads "${merged_filtered_rep_seqs}" \
            --p-n-jobs "${nthreads_classify}" \
            --o-classification "${taxonomy}"
    fi
    ### 9. Align filtered representative sequences (filtered-rep-seqs) to a tree.
    echo "$(date +"%F %H:%M:%S")"
    if [[ ! -f "${alignment}" ]]; then
        echo "Aligning the filtered representative sequences to a tree"
        qiime phylogeny align-to-tree-mafft-fasttree \
            --p-n-threads "${nthreads_classify}" \
            --i-sequences "${merged_filtered_rep_seqs}" \
            --o-alignment "${alignment}" \
            --o-masked-alignment "${masked_alignment}" \
            --o-tree "${unrooted_tree}" \
            --o-rooted-tree "${rooted_tree}"
    fi

done

echo "$(date +"%F %H:%M:%S")"