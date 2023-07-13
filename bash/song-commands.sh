qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filenames-song.tsv \
  --output-path song-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 

qiime demux summarize \
  --i-data song-paired-end-demux.qza \
  --o-visualization song-paired-end-demux.qzv
  
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs song-paired-end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table song-table-dada2-no-trunc.qza \
  --o-representative-sequences song-rep-seqs-dada2-no-trunc.qza \
  --o-denoising-stats song-stats-dada2-no-trunc.qza \
  --p-n-threads 0

######
# visualize your dada2 stats
qiime metadata tabulate \
  --m-input-file song-stats-dada2-no-trunc.qza \
  --o-visualization song-stats-dada2-no-trunc.qzv

  
# generate feature table (frequencies)
qiime feature-table summarize \
  --i-table song-table-dada2-no-trunc.qza \
  --m-sample-metadata-file filenames-song.tsv \
  --o-visualization song-table-dada2-no-trunc.qzv


# generate feature table (sequences)
qiime feature-table tabulate-seqs \
  --i-data song-rep-seqs-dada2-no-trunc.qza \
  --o-visualization song-rep-seqs-dada2-no-trunc.qzv 


###### Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences song-rep-seqs-dada2-no-trunc.qza \
  --o-alignment song-aligned-rep-seqs-dada2-no-trunc.qza \
  --o-masked-alignment song-masked-aligned-rep-seqs-dada2-no-trunc.qza \
  --o-tree song-unrooted-tree-dada2-no-trunc.qza \
  --o-rooted-tree song-rooted-tree-dada2-no-trunc.qza
  
  

##### RESCRIPT easy way	
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
  
