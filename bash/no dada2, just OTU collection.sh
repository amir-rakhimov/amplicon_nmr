# no dada2, just OTU collection
# dereplicating sequences
qiime vsearch dereplicate-sequences \
  --i-sequences single-end-demux.qza \
  --o-dereplicated-table table-vsearch.qza \
  --o-dereplicated-sequences rep-seqs-vsearch.qza
  
# de novo clustering
qiime vsearch cluster-features-de-novo \
  --i-table table-vsearch.qza \
  --i-sequences rep-seqs-vsearch.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-vsearch-dn-97.qza \
  --o-clustered-sequences rep-seqs-vsearch-dn-97.qza
  
# chimera checking
qiime vsearch uchime-denovo \
  --i-table table-vsearch-dn-97.qza \
  --i-sequences rep-seqs-vsearch-dn-97.qza \
  --output-dir uchime-vsearch-dn-out-97

# visualize results
qiime metadata tabulate \
  --m-input-file uchime-vsearch-dn-out-97/stats.qza \
  --o-visualization uchime-vsearch-dn-out-97/stats.qzv

# exclude chimeras but retain borderline chimeras
qiime feature-table filter-features \
  --i-table table-vsearch-dn-97.qza \
  --m-metadata-file uchime-vsearch-dn-out-97/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-table uchime-vsearch-dn-out-97/table-nonchimeric-w-borderline.qza
qiime feature-table filter-seqs \
  --i-data rep-seqs-vsearch-dn-97.qza \
  --m-metadata-file uchime-vsearch-dn-out-97/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-data uchime-vsearch-dn-out-97/rep-seqs-nonchimeric-w-borderline.qza
qiime feature-table summarize \
  --i-table uchime-vsearch-dn-out-97/table-nonchimeric-w-borderline.qza \
  --o-visualization uchime-vsearch-dn-out-97/table-nonchimeric-w-borderline.qzv