#!/usr/bin/env bash
# Obtain metadata, SRA files, and FASTQ files for each project.
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
sratoolkit_scripts_dir=~/sratoolkit.3.2.0-ubuntu64/bin
project_home_dir=~/projects/amplicon_nmr
bensch_metadata_folder="${project_home_dir}"/data/metadata/bensch-metadata
bensch_sra_dir="${project_home_dir}"/data/sra-files/bensch-sra
bensch_fastq_raw_dir="${project_home_dir}"/data/fastq/bensch-fastq/raw
##
liu_metadata_folder="${project_home_dir}"/data/metadata/liu-metadata
liu_sra_dir="${project_home_dir}"/data/sra-files/liu-sra
liu_fastq_raw_dir="${project_home_dir}"/data/fastq/liu-fastq/raw
##
sibai_metadata_folder="${project_home_dir}"/data/metadata/sibai-metadata
sibai_sra_dir="${project_home_dir}"/data/sra-files/sibai-sra
sibai_fastq_raw_dir="${project_home_dir}"/data/fastq/sibai-fastq/raw
##
shanmuganandam_metadata_folder="${project_home_dir}"/data/metadata/shanmuganandam-metadata
shanmuganandam_sra_dir="${project_home_dir}"/data/sra-files/shanmuganandam-sra
shanmuganandam_fastq_raw_dir="${project_home_dir}"/data/fastq/shanmuganandam-fastq/raw
mkdir -p "${project_home_dir}"/data/metadata/bensch-metadata
mkdir -p "${project_home_dir}"/data/sra-files/bensch-sra
mkdir -p "${project_home_dir}"/data/fastq/bensch-fastq/raw
##
mkdir -p "${project_home_dir}"/data/metadata/liu-metadata
mkdir -p "${project_home_dir}"/data/sra-files/liu-sra
mkdir -p "${project_home_dir}"/data/fastq/liu-fastq/raw

conda activate qc-tools
# Install SRA-Toolkit
curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    --output ~/sratoolkit-current.tar.gz
gunzip ~/sratoolkit-current.tar.gz
tar -xvf ~/sratoolkit-current.tar -C ~/sratoolkit
mv ~/sratoolkit/* ~/
rmdir ~/sratoolkit/

# Bensch data (Damaraland mole-rats)
# 1. Download SRA metadata (SRR IDs)
curl https://raw.githubusercontent.com/HannaBensch/FreezeDriedVSFrozen/refs/heads/main/data/FDvsFrozenMetadata.csv \
  --output "${bensch_metadata_folder}"/bensch-frozen-vs-freeze-dried.csv

esearch -db sra -query PRJNA781121 |efetch -format runinfo > "${bensch_metadata_folder}"/bensch-sra-runinfo.csv

# 2. Extract SRA Runs and BioSample columns into a separate file. Output is a comma-separated file.
# Run is column #1 (SRR16997153)
# BioSample is column #26 (SAMN23246763)
awk -F',' 'BEGIN{OFS=",";} NR==1 {print "SRA,BioSample"} NR>1 {print $1,$26}' \
  "${bensch_metadata_folder}"/bensch-sra-runinfo.csv > \
  "${bensch_metadata_folder}"/bensch-sra_run-biosample.csv 

# 3. Join SRA Runs and metadata into one file. Final columns are:
# Plate_No,sample,SampleNumber,NewSampleNumber,SampleOrder,Treatment,SampleDate,BioSample,Run
# 1,FD200,200,1,1,Freeze-dried,2019-04-01,SAMN23246723,SRR16997113

# Method: SRA and BioSample file will be a 'lookup' file
# sample_array[$2] is the common column (BioSample).
# Key is the BioSample, value is the SRA Run.
# If the 8-th column of the metadata file (i.e. BioSample) is in the array key, it'll print
# all of the metadata columns and the value of the array that corresponds to the
#  8-th column (sample_array[$8], i.e. SRA Run).
awk 'BEGIN {FS=OFS=","}
        NR==FNR {sample_array[$2]=$1;next} $8 in sample_array {print $0,sample_array[$8]}' \
    "${bensch_metadata_folder}"/bensch-sra_run-biosample.csv \
    "${bensch_metadata_folder}"/bensch-frozen-vs-freeze-dried.csv > \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-combined.csv

# 4. Filter by column #6 to keep only Frozen samples: awk '$COLUMN ~ /PATTERN/' file
awk 'BEGIN{FS=OFS=",";}
        NR ==1 ||$6 ~ /Frozen/' "${bensch_metadata_folder}"/bensch-sra_run-metadata-combined.csv > \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv

# 5. Extract sample names for the final metadata file (tab-separated). Final columns:
# sample-id	class	animal	sex	birthday	forward-absolute-filepath	reverse-absolute-filepath
# F200	DMR	Fukomys Damarensis	-	-	
# /home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F200_R1.fastq.gz	
# /home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F200_R2.fastq.gz

# Quote Fukomys Damarensis by adding  \" around the text
awk -F',' 'BEGIN{OFS="\t";} 
        NR==1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "sample-id","class","animal", "sex","birthday",
                    "forward-absolute-filepath", "reverse-absolute-filepath";
        print "";} 
        NR>1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "F"$3, "DMR","\"Fukomys Damarensis\"", "-", "-",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R1.fastq.gz",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R2.fastq.gz";
        print "";}' \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv > \
    "${bensch_metadata_folder}"/bensch-metadata.tsv

# 6. Run sratoolkit to download SRA files
# awk opens the file with metadata (only frozen samples) and extracts the column
# with SRA run names. The output is piped to prefetch with xargs -I {}
awk -F',' 'NR>1{print $9}' "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv | \
    xargs -I {}  "${sratoolkit_scripts_dir}"/prefetch {} --output-directory "${bensch_sra_dir}"
# "${sratoolkit_scripts_dir}"/prefetch  --option-file sra-list.txt --output-directory "${bensch_sra_dir}"

# For fastq-dump, take the metadata file, extract column 9 which contains SRA IDs, and 
# append the directory where you want to save them. Store awk output in an Array
mapfile -t sra_files_arr < <(awk -v sra_dir="${bensch_sra_dir}" -F',' \
    'NR>1{print sra_dir "/" $9 "/" $9 ".sra"}' \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv)
# mapfile -t files < <(awk -v dir="$DIR" '{print dir "/" $1 "/" $1 ".sra"}' file_list.txt)
# Uses awk to construct file paths.
# mapfile -t files stores output line by line into the array files.
# fastq-dump --split-files "${sra_files_arr[@]}"
# Expands the array to pass all files as arguments to fastq-dump.

# 7. Now use the array with fastq-dump
"${sratoolkit_scripts_dir}"/fastq-dump -v --outdir "${bensch_fastq_raw_dir}" --gzip --skip-technical \
    --readids --read-filter pass \
    --dumpbase --split-3 --clip "${sra_files_arr[@]}" 2>&1 | \
    tee "${project_home_dir}"/"${date_time}"_bensch_fastq_dump_report.txt


# sample-id	class	animal	sex	birthday	absolute-filepath
# #q2:types	categorical	categorical	categorical	categorical	categorical
# 2D10	NMR	naked mole rat	M	2016-11-23	/home/rakhimov/yasuda/data/raw/2D10_R1.fastq.gz

# 8. Rename SRA files to sample names
awk  -v fastq_dir="${bensch_fastq_raw_dir}" -F',' 'BEGIN{OFS="\n";} 
        NR>1 {printf  "%s\t%s",fastq_dir "/" $9 "_pass_1.fastq.gz",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R1.fastq.gz";
            print "";
            printf  "%s\t%s",fastq_dir "/" $9 "_pass_2.fastq.gz",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R2.fastq.gz";
            print "";}' \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv |\
    
    while read -r old new; do
    [[ -e "$old" ]] && mv "$old" "$new" && echo "Renamed: $old -> $new" || echo "File not found: $old"
done

    
# Liu
# 1. Download SRA Run information
esearch -db sra -query PRJNA428269 |efetch -format runinfo > "${liu_metadata_folder}"/liu-sra-runinfo.csv

# 2. Extract SRA Runs, sample names, and BioSample columns into a separate file. Output is a comma-separated file.
awk -F',' 'BEGIN{OFS=",";} NR==1 {print "SRA,BioSample,SampleName"} NR>1 {print $1,$26,$12}' \
  "${liu_metadata_folder}"/liu-sra-runinfo.csv > \
  "${liu_metadata_folder}"/liu-sra_run-biosample.csv 
# Run is column #1 (SRR6837860)
# BioSample is column #26 (SAMN08286762)
# Sample is column 12 (PAL_3)

# 3. Filter by #3 to keep only PVO samples: awk '$COLUMN ~ /PATTERN/' file
awk 'BEGIN{FS=OFS=",";}
        NR ==1 ||$3 ~ /PVO/' "${liu_metadata_folder}"/liu-sra_run-biosample.csv > \
    "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-only.csv 

# Add zero before single digits in sample names (e.g., PVO_1 -> PVO_01)
perl -pe 's/PVO_([1-9])\b/PVO_0\1/g' \
     "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-only.csv > \
     "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-digits-fix.csv
# perl → Runs the Perl interpreter.
# -p → Reads each line, processes it, and prints the modified version automatically.
# -e → Allows execution of the script provided as a command ('...').
# input.txt → The file being modified.
# s/.../.../g (Substitution Command)
    # s/// → A Perl substitution regex (s stands for substitute).
    # Left side (PVO_([1-9])\b) → The pattern to find.
    # Right side (PVO_0\1) → The replacement.
# Breaking Down the Pattern: PVO_([1-9])\b
    # PVO_ →  Matches the literal string "PVO_".
    # ([1-9]) →  Captures a single digit (1–9) inside ().
    # \1 in the replacement refers to this captured number.
    # \b (Word Boundary)
    # Ensures the match only occurs at the end of the word.
    # Prevents matching "PVO_10", "PVO_12", etc.
    # \b matches:
    # A space (PVO_1 → ✅ match)
    # A comma or period (PVO_1, → ✅ match)
    # End of line (PVO_1\n → ✅ match)
    # But NOT PVO_10 (no boundary before 0
# Breaking Down the Replacement: PVO_0\1
#     PVO_0 → Adds a leading zero.
#     \1 → Inserts the captured digit (1-9).
# g Flag (s/PVO_([1-9])\b/PVO_0\1/g)
# g → Global replacement (replaces all occurrences in a line, not just the first).

# 4. Filter out samples that have <20000 reads (identified beforehand)
awk 'BEGIN{FS=OFS=",";}
        NR ==1 ||$3 != "PVO_01"&& $3 != "PVO_02"&& $3 !=  "PVO_03" && $3 != "PVO_04" && \
            $3 != "PVO_05"&& $3 != "PVO_06"&& $3 != "PVO_11"&& $3 != "PVO_12"&& $3 != "PVO_13"' \
            "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-digits-fix.csv > \
    "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-for-qiime2.csv 

# Sort the file by sample name
(head -n 1 "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-for-qiime2.csv && \
    tail -n +2 "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-for-qiime2.csv | sort -t, -k3,3) > \
    "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-for-qiime2-sorted.csv

# 5. Extract sample names for the final metadata file (tab-separated). Final columns:
# sample-id  class  animal  sex  birthday  forward-absolute-filepath  reverse-absolute-filepath
# PVO_19  PVO  "Pteromys volans orii" -  -  
# /home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/PVO_19_R1.fastq.gz
# /home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/PVO_19_R2.fastq.gz

# Quote Fukomys Damarensis by adding  \" around the text
awk -F',' 'BEGIN{OFS="\t";} 
        NR==1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "sample-id","class","animal", "sex","birthday",
                    "forward-absolute-filepath", "reverse-absolute-filepath";
        print "";} 
        NR>1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", $3, "PVO","\"Pteromys volans orii\"", "-", "-",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/"$3"_R1.fastq.gz",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/"$3"_R2.fastq.gz";
        print "";}' \
    "${liu_metadata_folder}"/liu-sra_run-biosample-pvo-for-qiime2-sorted.csv > \
    "${liu_metadata_folder}"/liu-metadata.tsv








TODO






# 6. Run sratoolkit to download SRA files
# awk opens the file with metadata (only frozen samples) and extracts the column
# with SRA run names. The output is piped to prefetch with xargs -I {}
awk -F',' 'NR>1{print $9}' "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv | \
    xargs -I {}  "${sratoolkit_scripts_dir}"/prefetch {} --output-directory "${bensch_sra_dir}"
# "${sratoolkit_scripts_dir}"/prefetch  --option-file sra-list.txt --output-directory "${bensch_sra_dir}"

# For fastq-dump, take the metadata file, extract column 9 which contains SRA IDs, and 
# append the directory where you want to save them. Store awk output in an Array
mapfile -t sra_files_arr < <(awk -v sra_dir="${bensch_sra_dir}" -F',' \
    'NR>1{print sra_dir "/" $9 "/" $9 ".sra"}' \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv)
# mapfile -t files < <(awk -v dir="$DIR" '{print dir "/" $1 "/" $1 ".sra"}' file_list.txt)
# Uses awk to construct file paths.
# mapfile -t files stores output line by line into the array files.
# fastq-dump --split-files "${sra_files_arr[@]}"
# Expands the array to pass all files as arguments to fastq-dump.

# 7. Now use the array with fastq-dump
"${sratoolkit_scripts_dir}"/fastq-dump -v --outdir "${bensch_fastq_raw_dir}" --gzip --skip-technical \
    --readids --read-filter pass \
    --dumpbase --split-3 --clip "${sra_files_arr[@]}" 2>&1 | \
    tee "${project_home_dir}"/"${date_time}"_bensch_fastq_dump_report.txt


# sample-id	class	animal	sex	birthday	absolute-filepath
# #q2:types	categorical	categorical	categorical	categorical	categorical
# 2D10	NMR	naked mole rat	M	2016-11-23	/home/rakhimov/yasuda/data/raw/2D10_R1.fastq.gz

# 8. Rename SRA files to sample names
awk  -v fastq_dir="${bensch_fastq_raw_dir}" -F',' 'BEGIN{OFS="\n";} 
        NR>1 {printf  "%s\t%s",fastq_dir "/" $9 "_pass_1.fastq.gz",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R1.fastq.gz";
            print "";
            printf  "%s\t%s",fastq_dir "/" $9 "_pass_2.fastq.gz",
                    "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R2.fastq.gz";
            print "";}' \
    "${bensch_metadata_folder}"/bensch-sra_run-metadata-frozen-only.csv |\
    
    while read -r old new; do
    [[ -e "$old" ]] && mv "$old" "$new" && echo "Renamed: $old -> $new" || echo "File not found: $old"
done

# Sibai

# Shanmuganandam
