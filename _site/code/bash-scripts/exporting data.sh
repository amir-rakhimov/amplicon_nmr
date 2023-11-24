#do it  on your pc after you ran qiime2
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

# filter based on classes: keep only custom classes
qiime feature-table filter-samples \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}.qza \
  --m-metadata-file ${metadata_file} \
  --p-where "$query" \
  --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza
# visualisation: table of samples and features
qiime feature-table summarize \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza \
  --m-sample-metadata-file ${metadata_file} \
  --o-visualization pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qzv


# remove specific samples
echo SampleID > samples-to-remove.tsv
echo -e PVO_01\\nPVO_02\\nPVO_03\\nPVO_04\\nPVO_05\\nPVO_06\\nPVO_11\\nPVO_12\\nPVO_13>> samples-to-remove.tsv
# exclude-ids means we remove samples based on m-metadata-file
qiime feature-table filter-samples \
  --i-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza\
  --m-metadata-file samples-to-remove.tsv \
  --p-exclude-ids True \
  --o-filtered-table pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}-less_samples.qza
# visualisation of feature table
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





#if no filtering, use pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}.qza
# and pooled-single-filtered-rep-seqs-trimmed-dada2-234-finalfiltercheck.qzv
cp pooled-single-filtered-table-trimmed-dada2-234.qza pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza 

qiime tools export \
  --input-path pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.qza \
  --output-path  pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}

cd pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}/
mv feature-table.biom pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.biom  
mv pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.biom ../
cd ../
rm -rf pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}/

# we need to open representative sequences with qiime2view and save them as fasta file
# if sequences weren't filtered by sample, use pooled-single-filtered-rep-seqs-trimmed-dada2-234-finalfiltercheck.qzv
mv sequences.fasta pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.fasta


scp pooled-single-filtered-table-trimmed-dada2-${trunc_len}-${custom_classes_filename}.biom rakhimov@gw.ddbj.nig.ac.jp:~/picrust
scp pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-${custom_classes_filename}.fasta rakhimov@gw.ddbj.nig.ac.jp:~/picrust





../../../../picrust2-2.5.2/scripts/add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz                                             
../../../../picrust2-2.5.2/scripts/add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz                                             
../../../../picrust2-2.5.2/scripts/add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o pathways_out/path_abun_unstrat_descrip.tsv.gz  