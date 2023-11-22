# RESCRIPT hard way
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path tax_slv_ssu_132.txt \
  --output-path taxranks-silva-132-ssu-nr99.qza

qiime tools import \
  --type 'FeatureData[TaxidMap]' \
  --input-path taxmap_slv_ssu_ref_nr_132.txt \
  --output-path taxmap-silva-132-ssu-nr99.qza
  
qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path tax_slv_ssu_132.tre \
    --output-path taxtree-silva-132-nr99.qza 
	
# RESCRIPT easy way	
qiime rescript get-silva-data \
    --p-version '132' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-132-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-132-ssu-nr99-tax.qza

# RNA to DNA
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-132-ssu-nr99-rna-seqs.qza 
  --o-dna-sequences silva-132-ssu-nr99-seqs.qza
	
# remove sequences that contain 5 or more ambiguous bases (IUPAC compliant ambiguity bases) 
# and any homopolymers that are 8 or more bases in length. These are the default parameters. See the --help text for more details.

qiime rescript cull-seqs \
  --i-sequences silva-132-ssu-nr99-seqs.qza \
  --o-clean-sequences silva-132-ssu-nr99-seqs-cleaned.qza
	
# Extract reads
qiime feature-classifier extract-reads \
  --i-sequences silva-132-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
  --p-r-primer GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-132-ssu-nr99-seqs-cleaned-341f-805r.qza
  
# dereplicate extracted region  
qiime rescript dereplicate \
  --i-sequences silva-132-ssu-nr99-seqs-cleaned-341f-805r.qza \
  --i-taxa silva-132-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-132-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --o-dereplicated-taxa  silva-132-ssu-nr99-tax-341f-805r-derep-uniq.qza
	

# train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-132-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --i-reference-taxonomy silva-132-ssu-nr99-tax-341f-805r-derep-uniq.qza \
  --o-classifier silva-132-ssu-nr99-classifier-341f-805r-derep-uniq.qza
  
# test the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva-132-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads song-rep-seqs-dada2-no-trunc.qza \
  --o-classification song-taxonomy-dada2-no-trunc.qza

qiime metadata tabulate \
  --m-input-file song-taxonomy-dada2-no-trunc.qza \
  --o-visualization song-taxonomy-dada2-no-trunc.qzv

# taxa bar plots
qiime taxa barplot \
  --i-table song-table-dada2-no-trunc.qza \
  --i-taxonomy song-taxonomy-dada2-no-trunc.qza \
  --m-metadata-file filenames-song.tsv \
  --o-visualization song-taxa-bar-plots-dada2-no-trunc.qzv