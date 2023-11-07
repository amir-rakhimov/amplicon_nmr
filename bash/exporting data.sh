trunc_len=234
metadata_file="../../filenames-single-pooled-raw-supercomp.tsv"

custom_classes=("NMR" "B6mouse" "MSMmouse" "FVBNmouse" "DMR" "hare" "rabbit" "spalax"  "pvo")
custom_classes_as_string=$(IFS=,; echo "${custom_classes[*]}")
custom_classes_filename=$(IFS=-; echo "${custom_classes[*]}")

query=""
for item in "${custom_classes[@]}"; do
    query+="([class]='$item') OR "
done
query="${query% OR *}"

# filter based on classes
qiime feature-table filter-samples \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --m-metadata-file ${metadata_file} \
  --p-where "$query" \
  --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza

qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qzv


# remove specific samples
echo SampleID > samples-to-remove.tsv
echo -e PVO_01\\nPVO_02\\nPVO_03\\nPVO_04\\nPVO_05\\nPVO_06\\nPVO_11\\nPVO_12\\nPVO_13>> samples-to-remove.tsv

qiime feature-table filter-samples \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza\
  --m-metadata-file samples-to-remove.tsv \
  --p-exclude-ids True \
  --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}-less_samples.qza

qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}-less_samples.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}-less_samples.qzv
  
mv pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}-less_samples.qza  pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza


# filter representative sequences based on metadata(here we use feature table)
qiime feature-table filter-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza \
  --p-where "$query" \
  --p-no-exclude-ids \
  --o-filtered-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza
  
qiime feature-table tabulate-seqs \
  --i-data pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza \
  --o-visualization pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qzv





#if no filtering, use pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-finalfiltercheck.qzv


qiime tools export \
  --input-path pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza \
  --output-path  pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}

cd pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}/
mv feature-table.biom pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.biom  
mv pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.biom ../
cd ../
rm -rf pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}/


mv sequences.fasta pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.fasta


scp pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.biom rakhimov@gw.ddbj.nig.ac.jp:~/picrust
scp pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.fasta rakhimov@gw.ddbj.nig.ac.jp:~/picrust





../../../../picrust2-2.5.2/scripts/add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz                                             
../../../../picrust2-2.5.2/scripts/add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz                                             
../../../../picrust2-2.5.2/scripts/add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o pathways_out/path_abun_unstrat_descrip.tsv.gz  