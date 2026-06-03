#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.2
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "$(date +"%F %H:%M:%S")"

### 5. Create the reference database with RESCRIPT for taxonomic classification.
### Use SILVA v138.1
echo "$(date +"%F %H:%M:%S")"
if [[ -f "${silva_db_dir}"/rna-seqs.qza ]];then
    echo "${silva_db_dir}"/rna-seqs.qza exists
else
    echo "Creating the reference database with RESCRIPT for taxonomic classification"
    qiime rescript get-silva-data \
        --p-version "${silva_ver}" \
        --p-target 'SSURef_NR99' \
        --p-include-species-labels \
        --o-silva-sequences "${silva_db_dir}"/rna-seqs.qza \
        --o-silva-taxonomy "${silva_db_dir}"/tax.qza      
fi

echo "$(date +"%F %H:%M:%S")"
if [[ ! -f "${silva_db_dir}"/dna-seqs.qza ]];then
    echo "Reverse transcribing RNA to DNA"
    qiime rescript reverse-transcribe \
        --i-rna-sequences "${silva_db_dir}"/rna-seqs.qza \
        --o-dna-sequences "${silva_db_dir}"/dna-seqs.qza
fi

if [[ ! -f "${silva_db_dir}"/seqs-cleaned.qza ]]; then
    echo "Cleaning sequences: remove sequences that contain 5 or more ambiguous bases \
        (IUPAC compliant ambiguity bases) \
        and any homopolymers that are 8 or more bases in length. \
        These are the default parameters. See the --help text for more details."
    echo "$(date +"%F %H:%M:%S")"
    qiime rescript cull-seqs \
        --i-sequences "${silva_db_dir}"/dna-seqs.qza \
        --o-clean-sequences "${silva_db_dir}"/seqs-cleaned.qza
fi

echo "$(date +"%F %H:%M:%S")"
if [[ -f "${silva_db_dir}"/db-extracted-reads.qza ]];then
    echo "${silva_db_dir}"/db-extracted-reads.qza exists
else
    echo "Extracting reads"
    qiime feature-classifier extract-reads \
        --i-sequences "${silva_db_dir}"/seqs-cleaned.qza \
        --p-f-primer "${fwd_primer}" \
        --p-r-primer "${rev_primer}" \
        --p-trunc-len "${ref_db_trunc_len}" \
        --p-n-jobs "${nthreads_train_classifier}" \
        --p-read-orientation 'forward' \
        --o-reads "${silva_db_dir}"/db-extracted-reads.qza
fi

echo "$(date +"%F %H:%M:%S")"
if [[ -f "${silva_db_dir}"/tax-uniq.qza ]];then
    echo "${silva_db_dir}"/tax-uniq.qza exists
else
    echo "Dereplicating the extracted region"
    qiime rescript dereplicate \
        --i-sequences "${silva_db_dir}"/db-extracted-reads.qza \
        --i-taxa "${silva_db_dir}"/tax.qza \
        --p-mode 'uniq' \
        --o-dereplicated-sequences "${silva_db_dir}"/db-extracted-reads-uniq.qza \
        --o-dereplicated-taxa "${silva_db_dir}"/tax-uniq.qza
fi
### 6. Fit a Naive Bayes classifier
echo "$(date +"%F %H:%M:%S")"
if [[ -f "${silva_data_dir}"/naive-bayes-classifier.qza ]]; then
    echo "${silva_data_dir}"/naive-bayes-classifier.qza exists
else
    echo "Training a Naive Bayes classifier"
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "${silva_db_dir}"/db-extracted-reads-uniq.qza \
        --i-reference-taxonomy "${silva_db_dir}"/tax-uniq.qza \
        --o-classifier "${silva_data_dir}"/naive-bayes-classifier.qza
fi
echo "$(date +"%F %H:%M:%S")"