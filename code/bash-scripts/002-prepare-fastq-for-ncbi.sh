#!/usr/bin/env bash
# Rename NMR and B6 mouse 16S files
project_home_dir=~/projects/amplicon_nmr
nmr_b6_16s_original_dir="${project_home_dir}"/data/fastq/for-ncbi/"NMR B6 microbiome_16S_noster"/fastq
nmr_b6_16s_renamed_dir="${project_home_dir}"/data/fastq/for-ncbi/renamed

msm_fvbn_16s_original_dir="${project_home_dir}"/data/fastq/for-ncbi/"MSM and FVB Fecal Miseq FastQ"
msm_fvbn_16s_renamed_dir="${project_home_dir}"/data/fastq/for-ncbi/renamed


for FILE in "${nmr_b6_16s_original_dir}"/*.fastq.gz
do
    # echo "renaming ${FILE}"
    SAMPLE=$(basename "${FILE}"| sed "s/kumamoto-//" | sed "s/_S[0-9]\{2\}_L001//")
    SAMPLE=$(echo "${SAMPLE}" | sed "s/_001//"| sed "s/mf-/mf_/" | sed "s/R1\.fastq\.gz/16S_R1\.fastq\.gz/" | sed "s/R2\.fastq\.gz/16S_R2\.fastq\.gz/")
    echo $(basename "${FILE}")
    echo "${SAMPLE}"
    cp "${FILE}" "${nmr_b6_16s_renamed_dir}"/"${SAMPLE}";
done

mv "${nmr_b6_16s_renamed_dir}"/DGC6_16S_R1.fastq.gz "${nmr_b6_16s_renamed_dir}"/DCG6_16S_R1.fastq.gz 
mv "${nmr_b6_16s_renamed_dir}"/DGC6_16S_R2.fastq.gz "${nmr_b6_16s_renamed_dir}"/DCG6_16S_R2.fastq.gz 

# Rename MSM and FVBN 16S files
for FILE in "${msm_fvbn_16s_original_dir}"/*.fastq.gz
do
    SAMPLE=$(basename "${FILE}" | sed "s/_S[0-9]\{2\}_L001//")
    SAMPLE=$(echo "${SAMPLE}" | sed "s/_001//"| sed "s/R1\.fastq\.gz/16S_R1\.fastq\.gz/" | sed "s/R2\.fastq\.gz/16S_R2\.fastq\.gz/")
    SAMPLE=$(echo "${SAMPLE}" | sed "s/^3/MSM3/" |  sed "s/^5/FVBN5/" )
    echo "${SAMPLE}"
    cp "${FILE}" "${msm_fvbn_16s_renamed_dir}"/"${SAMPLE}"
done



# Rename WMS files
# Rename files
cd "${fastq_dir}"
for FILE in nmrF*.fq.gz; do mv "${FILE}" $(echo "${FILE}" | sed 's/nmrF_//'); done
for FILE in *DKDN230040401-1A_HC3LYDSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040401-1A_HC3LYDSX7/wms/'); done
for FILE in *DKDN230040401-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040401-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040402-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040402-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040410-1A_HC52YDSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040410-1A_HC52YDSX7/wms/'); done
for FILE in *DKDN230040409-1A_HC3LYDSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040409-1A_HC3LYDSX7/wms/'); done
for FILE in *DKDN230040409-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040409-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040411-1A_HC52YDSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040411-1A_HC52YDSX7/wms/'); done
for FILE in *DKDN230040407-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040407-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040408-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040408-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040406-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040406-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040403-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040403-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040404-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040404-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040405-1A_HC3LYDSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040405-1A_HC3LYDSX7/wms/'); done
for FILE in *DKDN230040405-1A_HC557DSX7*.fq.gz; 
  do mv "${FILE}" $(echo "${FILE}" | sed 's/DKDN230040405-1A_HC557DSX7/wms/'); done
