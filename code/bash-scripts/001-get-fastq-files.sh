#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
# Obtain metadata, SRA files, and FASTQ files for each project.
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}

cd "${project_home_dir}"
for author in "${authors[@]}"; do
    echo "${author}"
    author_metadata_folder="${all_metadata_dir}"/"${author}"-metadata
    author_sra_dir="${sra_files_dir}"/"${author}"-sra
    author_fastq_raw_dir="${fastq_dir}"/"${author}"-fastq/raw
    author_sra_to_sample_map_file="${author_metadata_folder}"/"${author}"-sra_run-biosample-filtered.csv

    mkdir -p "${author_metadata_folder}"
    mkdir -p "${author_sra_dir}"
    mkdir -p "${author_fastq_raw_dir}"
    if [[ "${author}" == "bensch" ]]; then
        # 0. Download custom SRA metadata (SRR IDs)
        curl https://raw.githubusercontent.com/HannaBensch/FreezeDriedVSFrozen/refs/heads/main/data/FDvsFrozenMetadata.csv \
            --output "${author_metadata_folder}"/"${author}"-frozen-vs-freeze-dried.csv
        bioproject_id=PRJNA781121
    elif [[ "${author}" == "liu" ]]; then
        bioproject_id=PRJNA428269
    elif [[ "${author}" == "sibai" ]]; then
        bioproject_id=PRJNA607251
    elif [[ "${author}" == "shanmuganandam" ]]; then
        bioproject_id=PRJNA576096
    else
        sleep 1
    fi

    # 1. Download SRA Run information
    esearch -db sra -query "${bioproject_id}" |efetch -format runinfo > \
        "${author_metadata_folder}"/"${author}"-sra-runinfo.csv

    # 2. Extract SRA Runs, sample names, and BioSample columns into a separate file. 
    # Output is a comma-separated file.
    if [[ "${author}" == "shanmuganandam" ]]; then
        # 2. Extract SRA Runs, sample names, and BioSample columns into a separate file. Output is a comma-separated file.
        awk -F',' 'BEGIN{OFS=",";} NR==1  {print "SRA,BioSample,SampleName"} NR>1  && $0 ~ /Illumina/ {print $1,$26,$12}' \
        "${author_metadata_folder}"/"${author}"-sra-runinfo.csv \
        > "${author_metadata_folder}"/"${author}"-sra_run-biosample.csv 
    else    
        awk -F',' 'BEGIN{OFS=",";} NR==1 {print "SRA,BioSample,SampleName"} NR>1 {print $1,$26,$12}' \
        "${author_metadata_folder}"/"${author}"-sra-runinfo.csv > \
        "${author_metadata_folder}"/"${author}"-sra_run-biosample.csv 
        # Run is column #1 (SRR6837860)
        # BioSample is column #26 (SAMN08286762)
        # Sample is column 12 (522_library)
    fi

    if [[ "${author}" == "bensch" ]]; then
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
            "${author_metadata_folder}"/"${author}"-sra_run-biosample.csv \
            "${author_metadata_folder}"/"${author}"-frozen-vs-freeze-dried.csv | \
            # Filter to keep only frozen samples
            awk 'BEGIN{FS=OFS=",";}
                NR ==1 {print $9, $8, $3};
                NR>1 && $6 ~ /Frozen/ {print $9, $8, "F"$3}' |\
            # 4. Sort the file by sample name
            # The key is in the 1q -- print first line (header) and quit (leaving the rest of the input to sort). https://stackoverflow.com/questions/14562423/is-there-a-way-to-ignore-header-lines-in-a-unix-sort
            { sed -u 1q; sort -t, -k3,3; } \
                > "${author_sra_to_sample_map_file}"
        
        # 5. Extract sample names for the final metadata file (tab-separated). Final columns:
        # sample-id	class	animal	sex	birthday	forward-absolute-filepath	reverse-absolute-filepath
        # F200	DMR	Fukomys Damarensis	-	-	
        # /home/rakhimov/projects/amplicon_nmr/data/fastq/"${author}"-fastq/raw/F200_R1.fastq.gz	
        # /home/rakhimov/projects/amplicon_nmr/data/fastq/"${author}"-fastq/raw/F200_R2.fastq.gz

        # Quote Fukomys Damarensis by adding  \" around the text
        awk -F',' 'BEGIN{OFS="\t";} 
                NR==1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "sample-id","class","animal", "sex","birthday",
                            "forward-absolute-filepath", "reverse-absolute-filepath";
                print "";} 
                NR>1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "F"$3, "DMR","\"Fukomys Damarensis\"", "-", "-",
                            "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R1.fastq.gz",
                            "/home/rakhimov/projects/amplicon_nmr/data/fastq/bensch-fastq/raw/F"$3"_R2.fastq.gz";
                print "";}' \
            "${author_sra_to_sample_map_file}" > \
            "${author_metadata_folder}"/"${author}"-metadata.tsv
            
    elif [[ "${author}" == "liu" ]]; then
        # 3. Filter by #3 to keep only PVO samples: awk '$COLUMN ~ /PATTERN/' file
        awk 'BEGIN{FS=OFS=",";}
            NR ==1 ||$3 ~ /PVO/' "${author_metadata_folder}"/"${author}"-sra_run-biosample.csv | \
            # Add zero before single digits in sample names (e.g., S1 -> S01)
            perl -pe 's/PVO_([1-9])\b/PVO_0\1/g' | \
            # 4. Sort the file by sample name
            # The key is in the 1q -- print first line (header) and quit (leaving the rest of the input to sort). https://stackoverflow.com/questions/14562423/is-there-a-way-to-ignore-header-lines-in-a-unix-sort
            { sed -u 1q; sort -t, -k3,3; } |
            # Filter out samples that have <20000 reads (identified beforehand)
            awk 'BEGIN{FS=OFS=",";}
                NR ==1 ||$3 != "PVO_01"&& $3 != "PVO_02"&& $3 !=  "PVO_03" && $3 != "PVO_04" && \
                    $3 != "PVO_05"&& $3 != "PVO_06"&& $3 != "PVO_11"&& $3 != "PVO_12"&& $3 != "PVO_13"' \
            >  "${author_sra_to_sample_map_file}"
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
                
        # 5. Extract sample names for the final metadata file (tab-separated). Final columns:
        # sample-id  class  animal  sex  birthday  forward-absolute-filepath  reverse-absolute-filepath
        # PVO_19  PVO  "Pteromys volans orii" -  -  
        # /home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/PVO_19_R1.fastq.gz
        # /home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/PVO_19_R2.fastq.gz

        # Quote Pteromys volans orii by adding  \" around the text
        awk -F',' 'BEGIN{OFS="\t";} 
                NR==1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "sample-id","class","animal", "sex","birthday",
                            "forward-absolute-filepath", "reverse-absolute-filepath";
                print "";} 
                NR>1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", $3, "PVO","\"Pteromys volans orii\"", "-", "-",
                            "/home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/"$3"_R1.fastq.gz",
                            "/home/rakhimov/projects/amplicon_nmr/data/fastq/liu-fastq/raw/"$3"_R2.fastq.gz";
                print "";}' \
            "${author_sra_to_sample_map_file}" > \
            "${author_metadata_folder}"/"${author}"-metadata.tsv
    elif [[ "${author}" == "sibai" ]]; then
        # 3. Filter by #3 to keep only T0SpalaxFecal samples: awk '$COLUMN ~ /PATTERN/' file
        awk 'BEGIN{FS=OFS=",";}
                NR ==1 ||$3 ~ /T0SpalaxFecal/' "${author_metadata_folder}"/"${author}"-sra_run-biosample.csv |\
                sed 's/T0SpalaxFecal[0-9]*_//' | \
            # Add zero before single digits in sample names (e.g., S1 -> S01)
            perl -pe 's/S([1-9])\b/S0\1/g' | \
            # 4. Sort the file by sample name
            # The key is in the 1q -- print first line (header) and quit (leaving the rest of the input to sort). https://stackoverflow.com/questions/14562423/is-there-a-way-to-ignore-header-lines-in-a-unix-sort
            { sed -u 1q; sort -t, -k3,3; } \
            >  "${author_sra_to_sample_map_file}"
        
        # 5. Extract sample names for the final metadata file (tab-separated). Final columns:
        # sample-id  class  animal  sex  birthday  forward-absolute-filepath  reverse-absolute-filepath
        # S01  spalax  "Nannospalax leucodon" -  -  
        # /home/rakhimov/projects/amplicon_nmr/data/fastq/sibai-fastq/raw/S01_R1.fastq.gz
        # /home/rakhimov/projects/amplicon_nmr/data/fastq/sibai-fastq/raw/S01_R2.fastq.gz

        # Quote Nannospalax leucodon by adding  \" around the text
        species_name="Nannospalax leucodon"
        class_name="spalax"
        awk -v author="${author}" -v species_name="${species_name}" -v class_name="${class_name}" \
            -F',' 'BEGIN{OFS="\t";} 
                NR==1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", "sample-id","class","animal", "sex","birthday",
                            "forward-absolute-filepath", "reverse-absolute-filepath";
                print "";} 
                NR>1 {printf  "%s\t%s\t%s\t%s\t%s\t%s\t%s", $3, class_name,"\"Nannospalax leucodon\"", "-", "-",
                            "/home/rakhimov/projects/amplicon_nmr/data/fastq/" author "-fastq/raw/"$3"_R1.fastq.gz",
                            "/home/rakhimov/projects/amplicon_nmr/data/fastq/author-fastq/raw/"$3"_R2.fastq.gz";
                print "";}' \
             "${author_sra_to_sample_map_file}" > \
            "${author_metadata_folder}"/"${author}"-metadata.tsv
    elif [[ "${author}" == "shanmuganandam" ]]; then
        # 3. Filter by #3 to keep only T0SpalaxFecal samples: awk '$COLUMN ~ /PATTERN/' file
        awk 'BEGIN{FS=OFS=",";}
                NR ==1 || $3 == "MF_136" || $3 == "MF_139" || $3 == "MF_140" || $3 == "MF_143" || $3 == "MF_144" || 
                    $3 == "MF_147" || $3 == "MF_148" || $3 == "MF_149" || $3 == "MF_150" || $3 == "MF_149" || $3 == "MF_151" || 
                    $3 == "MF_152" || $3 == "MF_153" || $3 == "MF_154" || $3 == "MF_155" || $3 == "MF_156" ' \
                "${author_metadata_folder}"/"${author}"-sra_run-biosample.csv |\
            # 4. Sort the file by sample name
            # The key is in the 1q -- print first line (header) and quit (leaving the rest of the input to sort). https://stackoverflow.com/questions/14562423/is-there-a-way-to-ignore-header-lines-in-a-unix-sort
            { sed -u 1q; sort -t, -k3,3; } \
            >  "${author_sra_to_sample_map_file}"


    fi

    # Next steps are same for all samples

    # 6. Run sratoolkit to download SRA files
    # awk opens the file with metadata and extracts the column
    # with SRA run names. The output is piped to prefetch with xargs -I {}
    awk -F',' 'NR>1{print $1}' "${author_sra_to_sample_map_file}" | \
        xargs -I {}  "${sratoolkit_scripts_dir}"/prefetch {} --output-directory "${author_sra_dir}"  2>&1 | \
        tee "${project_home_dir}"/jobreports/"${date_time}"_"${author}"_prefetch_report.txt
    # "${sratoolkit_scripts_dir}"/prefetch  --option-file sra-list.txt --output-directory "${author_sra_dir}"

    # For fastq-dump, take the metadata file, extract column 1 which contains SRA IDs, and 
    # append the directory where you want to save them. Store awk output in an Array
    mapfile -t sra_files_arr < <(awk -v sra_dir="${author_sra_dir}" -F',' \
        'NR>1{print sra_dir "/" $1 "/" $1 ".sra"}' \
        "${author_sra_to_sample_map_file}")
    # mapfile -t files < <(awk -v dir="$DIR" '{print dir "/" $1 "/" $1 ".sra"}' file_list.txt)
    # Uses awk to construct file paths.
    # mapfile -t files stores output line by line into the array files.
    # fastq-dump --split-files "${sra_files_arr[@]}"
    # Expands the array to pass all files as arguments to fastq-dump.

    # 7. Now use the array with fastq-dump
    "${sratoolkit_scripts_dir}"/fastq-dump -v --outdir "${author_fastq_raw_dir}" --gzip --skip-technical \
        --readids --read-filter pass \
        --dumpbase --split-3 --clip "${sra_files_arr[@]}" 2>&1 | \
        tee "${project_home_dir}"/jobreports/"${date_time}"_"${author}"_fastq_dump_report.txt

    # sample-id	class	animal	sex	birthday	absolute-filepath
    # #q2:types	categorical	categorical	categorical	categorical	categorical
    # 2D10	NMR	naked mole rat	M	2016-11-23	/home/rakhimov/yasuda/data/raw/2D10_R1.fastq.gz

    # 8. Rename SRA files to sample names
    awk  -v fastq_dir="${author_fastq_raw_dir}" -v author="${author}" -F',' 'BEGIN{OFS="\n";} 
            NR>1 {printf  "%s\t%s",fastq_dir "/" $1 "_pass_1.fastq.gz",
                        "/home/rakhimov/projects/amplicon_nmr/data/fastq/" author "-fastq/raw/"$3"_R1.fastq.gz";
                print "";
                printf  "%s\t%s",fastq_dir "/" $1 "_pass_2.fastq.gz",
                        "/home/rakhimov/projects/amplicon_nmr/data/fastq/" author "-fastq/raw/"$3"_R2.fastq.gz";
                print "";}' \
        "${author_sra_to_sample_map_file}" |\
        
        while read -r old new; do
        [[ -e "$old" ]] && mv "$old" "$new" && echo "Renamed: $old -> $new" || echo "File not found: $old"
    done 2>&1 | \
        tee "${project_home_dir}"/jobreports/"${date_time}"_"${author}"_rename_report.txt



done
