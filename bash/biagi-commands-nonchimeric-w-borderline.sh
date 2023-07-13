#!/bin/bash
# import data
qiime tools import \
  --type  'SampleData[SequencesWithQuality]' \
  --input-path filenames.txt \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 

# quality scores: 1 722 650 reads in total
qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization single-end-demux.qzv
# Min: 48678; Max: 49499; Mean: 49218.571429

#####
# denoising with dada2
# attempt 1: trim before position 23, truncate after 455
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-demux.qza \
  --p-trim-left 23 \
  --p-trunc-len 455 \
  --o-representative-sequences rep-seqs-dada2-initial.qza \
  --o-table table-dada2-initial.qza \
  --o-denoising-stats stats-dada2-initial.qza
  --p-n-threads 8

# attempt 2: no truncation
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences rep-seqs-dada2-no-trunc.qza \
  --o-table table-dada2-no-trunc.qza \
  --o-denoising-stats stats-dada2-no-trunc.qza \
 # --p-n-threads 8  

# attempt 3: denoising with dada2 and truncating at 439
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 439 \
  --o-representative-sequences rep-seqs-dada2-439-trunc.qza \
  --o-table table-dada2-439-trunc.qza \
  --o-denoising-stats stats-dada2-439-trunc.qza \
#  --p-n-threads 8


######
# visualize your dada2 stats
qiime metadata tabulate \
  --m-input-file stats-dada2-no-trunc.qza \
  --o-visualization stats-dada2-no-trunc.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2-439-trunc.qza \
  --o-visualization stats-dada2-439-trunc.qzv
  
# generate feature table (frequencies)
qiime feature-table summarize \
  --i-table table-dada2-no-trunc.qza \
  --o-visualization table-dada2-no-trunc.qzv

qiime feature-table summarize \
  --i-table table-dada2-439-trunc.qza \
  --o-visualization table-dada2-439-trunc.qzv

# generate feature table (sequences)
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2-no-trunc.qza \
  --o-visualization rep-seqs-dada2-no-trunc.qzv 

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2-439-trunc.qza \
  --o-visualization rep-seqs-dada2-439-trunc.qzv 

#####
#clustering: not truncated
qiime vsearch cluster-features-de-novo \
  --i-table table-dada2-no-trunc.qza \
  --i-sequences rep-seqs-dada2-no-trunc.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-dada2-no-trunc-dn-97.qza \
  --o-clustered-sequences rep-seqs-dada2-no-trunc-dn-97.qza

## visualize clustering results
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2-no-trunc-dn-97.qza \
  --o-visualization rep-seqs-dada2-no-trunc-dn-97.qzv 

# chimera check: de novo with uchime
qiime vsearch uchime-denovo \
  --i-table table-dada2-no-trunc-dn-97.qza \
  --i-sequences rep-seqs-dada2-no-trunc-dn-97.qza \
  --output-dir uchime-dada2-no-trunc-dn-out-97
  
## visualize summary stats
qiime metadata tabulate \
  --m-input-file uchime-dada2-no-trunc-dn-out-97/stats.qza \
  --o-visualization uchime-dada2-no-trunc-dn-out-97/stats.qzv

# chimera filtering 1: retain borderline
qiime feature-table filter-features \
  --i-table table-dada2-no-trunc-dn-97.qza \
  --m-metadata-file uchime-dada2-no-trunc-dn-out-97/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-table uchime-dada2-no-trunc-dn-out-97/table-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline.qza
  
qiime feature-table filter-seqs \
  --i-data rep-seqs-dada2-no-trunc-dn-97.qza \
  --m-metadata-file uchime-dada2-no-trunc-dn-out-97/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-data uchime-dada2-no-trunc-dn-out-97/rep-seqs-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline.qza
  
qiime feature-table summarize \
  --i-table uchime-dada2-no-trunc-dn-out-97/table-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline.qza \
  --o-visualization uchime-dada2-no-trunc-dn-out-97/table-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline.qzv

#####
# taxonomic classification
# import otus
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ./gg_13_5_otus/rep_set/97_otus.fasta \
  --output-path gg_13_5_97_otus.qza

# import taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ./gg_13_5_otus/taxonomy/97_otu_taxonomy.txt \
  --output-path gg_13_5_97_ref-taxonomy.qza

# extract reference reads: no truncation, no min or max length filtering
qiime feature-classifier extract-reads \
  --i-sequences gg_13_5_97_otus.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len -1 \
  --p-min-length 0 \
  --p-max-length 0 \
  --o-reads gg_13_5_97_ref-seqs-no-minmax-no-len-filter.qza

# train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads gg_13_5_97_ref-seqs-no-minmax-no-len-filter.qza \
  --i-reference-taxonomy gg_13_5_97_ref-taxonomy.qza \
  --o-classifier gg_13_5_97_classifier-no-minmax-no-len-filter.qza

# test the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier gg_13_5_97_classifier-no-minmax-no-len-filter.qza \
  --i-reads uchime-dada2-no-trunc-dn-out-97/rep-seqs-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline.qza \
  --o-classification taxonomy-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline-no-minmax-no-len-filter.qza

# SILVA
qiime feature-classifier classify-sklearn \
  --i-classifier ../song/silva-132-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads rep-seqs-dada2-no-trunc.qza \
  --o-classification taxonomy-dada2-silva.qza


qiime metadata tabulate \
  --m-input-file taxonomy-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline-no-minmax-no-len-filter.qza \
  --o-visualization taxonomy-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline-no-minmax-no-len-filter.qzv

# taxa bar plots
qiime taxa barplot \
  --i-table uchime-dada2-no-trunc-dn-out-97/table-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline.qza \
  --i-taxonomy taxonomy-dada2-no-trunc-dn-97-uchime-nonchimeric-w-borderline-no-minmax-no-len-filter.qza \
  --m-metadata-file filenames.txt \
  --o-visualization taxa-bar-plots-dada2-no-trunc-dn-97-uchime-w-borderline-no-minmax-no-len-filter.qzv

#####
# generate a phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2-no-trunc.qza \
  --o-alignment aligned-rep-seqs-dada2-no-trunc.qza \
  --o-masked-alignment masked-aligned-rep-seqs-dada2-no-trunc.qza \
  --o-tree unrooted-tree-dada2-no-trunc.qza \
  --o-rooted-tree rooted-tree-dada2-no-trunc.qza

# alpha and beta diversity
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-dada2.qza \
  --m-metadata-file filenames.txt \
  --p-sampling-depth 3400 \
  --output-dir core-metrics-results

