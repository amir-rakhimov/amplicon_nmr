#!/usr/bin/env bash
conda activate qiime2
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138-1-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138-1-ssu-nr99-tax.qza
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138-1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138-1-ssu-nr99-seqs.qza
qiime rescript cull-seqs \
  --i-sequences silva-138-1-ssu-nr99-seqs.qza \
  --o-clean-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza
qiime feature-classifier extract-reads \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-n-jobs 2 \
  --p-read-orientation 'forward' \
  --o-reads silva-138-1-ssu-nr99-seqs-cleaned-341f-805r.qza
qiime rescript dereplicate \
  --i-sequences silva-138-1-ssu-nr99-seqs-cleaned-341f-805r.qza \
  --i-taxa silva-138-1-ssu-nr99-tax.qza \
  --p-rank-handles 'silva' \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-1-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --o-dereplicated-taxa  silva-138-1-ssu-nr99-tax-341f-805r-derep-uniq.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-1-ssu-nr99-seqs-cleaned-341f-805r-uniq.qza \
  --i-reference-taxonomy silva-138-1-ssu-nr99-tax-341f-805r-derep-uniq.qza \
  --o-classifier silva-138-1-ssu-nr99-classifier-341f-805r-derep-uniq.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-1-ssu-nr99-classifier-341f-805r-derep-uniq.qza \
  --i-reads yashliu-rep-seqs-trimmed-dada2-284-203.qza \ # not filtered \
  --o-classification yashliu-taxonomy-138-1-trimmed-dada2-284-203.qza
qiime metadata tabulate \
  --m-input-file yashliu-taxonomy-trimmed-dada2-284-203.qza \
  --o-visualization yashliu-taxonomy-138-1-trimmed-dada2-284-203.qzv
qiime taxa barplot \
  --i-table yashliu-table-trimmed-dada2-284-203.qza \
  --i-taxonomy yashliu-taxonomy-trimmed-dada2-284-203.qza \
  --m-metadata-file filenames-yashliu-raw-supercomp.tsv \
  --o-visualization yashliu-taxa-bar-plots-138-1-trimmed-dada2-284-203.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences yashliu-rep-seqs-trimmed-dada2-284-203.qza \
  --o-alignment yashliu-aligned-138-1-rep-seqs-trimmed-dada2-284-203.qza \
  --o-masked-alignment yashliu-masked-aligned-138-1-rep-seqs-trimmed-dada2-284-203.qza \
  --o-tree yashliu-unrooted-tree-138-1-trimmed-dada2-284-203.qza \
  --o-rooted-tree yashliu-rooted-tree-138-1-trimmed-dada2-284-203.qza