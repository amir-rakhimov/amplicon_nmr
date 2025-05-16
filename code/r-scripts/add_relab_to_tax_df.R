# The function add_relab_to_tax_df: 1. Takes a table with absolute abundances
# 2. Groups rows by class and samples
# 3. Then adds a column with total abundance (sum of all reads) per
# sample (TotalSample, which is a sum of Abundance column per sample per host)
# 4. Then groups by class, sample, and a the lowest taxonomic rank (tax.rank)
# 5. Then adds a column with relative abundances, which are absolute abundances
# of each taxon in a sample divided by the total number of reads in a sample,
# multiplied by 100.
# 6. Then, the dataframe is grouped by class 
# 7. And we add a TotalClass column, which is a sum of reads from all samples
# belonging to a class.
# 8. Then, we group by class and the lowest taxonomic rank (tax.rank).
# 9. We add a TotalAgglomRank column, which is a sum of all reads of each 
# tax.rank in each sample.
# 8. Finally, the MeanRelativeAbundance column is the total number of reads
# belonging to a certain taxon in a host divided by the total number of reads
# from that host, multiplied by 100.
# (MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
add_relab_to_tax_df<-function(tax.df,tax.rank){
  tax.df<-tax.df%>%
    group_by(class,Sample)%>%
    mutate(TotalSample=sum(Abundance))%>%
    group_by_at(c("class","Sample",tax.rank))%>%
    mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
    group_by(class)%>% # group by class (animal host),
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",tax.rank))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  return(tax.df)
}
