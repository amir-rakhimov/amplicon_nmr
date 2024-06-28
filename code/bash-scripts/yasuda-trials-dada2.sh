#!/usr/bin/env bash
trunc_len_f=313
trunc_len_r=229
front_primer=CCTACGGGNGGCWGCAG
rev_primer=GACTACHVGGGTATCTAATCC
front_primer_position=341f
rev_primer_position=805r
ref_db_trunc_len=465
#  --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
#   --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
trunc_len_f=284
trunc_len_r=229

trunc_len_f=284
trunc_len_r=203

trunc_len_f=273
trunc_len_r=203

trunc_len_f=273
trunc_len_r=198

trunc_len_f=265
trunc_len_r=203

trunc_len_f=265
trunc_len_r=198
conda activate qiime2
qiime tools import \
  --type  'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filenames-yasuda.tsv \
  --output-path yasuda-paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 
qiime demux summarize \
  --i-data yasuda-paired-end-demux.qza \
  --o-visualization yasuda-paired-end-demux.qzv
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences yasuda-paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f ${front_primer} \
  --p-front-r ${rev_primer} \
  --o-trimmed-sequences yasuda-trimmed-paired-end.qza
qiime demux summarize \
  --i-data yasuda-trimmed-paired-end.qza \
  --o-visualization yasuda-trimmed-paired-end.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs yasuda-trimmed-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f ${trunc_len_f} \
  --p-trunc-len-r ${trunc_len_r} \
  --o-table yasuda-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-representative-sequences yasuda-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-denoising-stats yasuda-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --p-n-threads 0
qiime metadata tabulate \
  --m-input-file yasuda-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization yasuda-stats-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-table summarize \
  --i-table yasuda-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization yasuda-table-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv
qiime feature-table tabulate-seqs \
  --i-data yasuda-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qza \
  --o-visualization yasuda-rep-seqs-trimmed-dada2-${trunc_len_f}-${trunc_len_r}.qzv 