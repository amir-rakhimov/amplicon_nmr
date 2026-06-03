#!/usr/bin/env bash
project_home_dir=$HOME/projects/amplicon_nmr
# read_type="single_end"
read_type=paired_end
### Truncation length for DADA2
# trunc_len_f=234

# trunc_len_f=284
# trunc_len_r=203
# trunc_len_f=290 # incorrect because we remove primers
# trunc_len_r=235 # incorrect because we remove primers
# trunc_len_f=273
# trunc_len_r=214
nthreads_qc=4
nthreads_dada2=16
nthreads_train_classifier=25
nthreads_classify=20
analysed_data=shared

# RUN_ID="$(date +%Y%m%d_%H%M%S)"-"${analysed_data}"-"${read_type}"
# RUN_ID=20260515_164224-shared-paired_end
# RUN_ID=20260518_214705-shared-paired_end
RUN_ID=20260521_191423-shared-paired_end
authors=("yasuda" "bensch" "liu" "shanmuganandam" "sibai")
final_authors=("yasuda" "bensch" "shanmuganandam" "sibai")
pairs=("278 229" "278 225" "278 220" "278 215" "275 229" "275 225" "270 229" "270 225" "265 220" "230 230" "230 225" "225 230" "225 225" "270 210" "259 199" "284 203" "273 214" ) 
final_pairs=("230 225" "225 225") 
# Output directory for this run
OUTDIR="${project_home_dir}"/results/1-qiime2/"${RUN_ID}"
# QIIME2 output directory
qiime2_output_dir="${OUTDIR}"/01-qiime2_output
# SILVA database files directory
silva_data_dir="${OUTDIR}"/02-silva_data
# FASTQ files directory
fastq_dir=data/0-raw/fastq
# QC output directory (only FastQC and MultiQC)
fastqc_raw_output_dir=data/1-qc/fastqc_raw/"${RUN_ID}"
multiqc_raw_output_dir=data/1-qc/multiqc_raw/"${RUN_ID}"
fastqc_trimmed_output_dir=data/1-qc/fastqc_trimmed/"${RUN_ID}"
multiqc_trimmed_output_dir=data/1-qc/multiqc_trimmed/"${RUN_ID}"
# Metadata file contains the fastq file paths
metadata_dir="${project_home_dir}"/data/metadata/"${analysed_data}"
metadata_file="${metadata_dir}"/filenames-paired.tsv

# SILVA database parameters
silva_ver="138.1" # for QIIME2 
silva_v_id=138_1 # for directory name
silva_db_dir="${silva_data_dir}"/"${silva_v_id}"-SSURef_NR99 # directory name

### Front and reverse primers
fwd_primer=CCTACGGGNGGCWGCAG
rev_primer=GACTACHVGGGTATCTAATCC

### Truncation length for the taxonomic reference database
ref_db_trunc_len=465
# ref_db_trunc_len=234

