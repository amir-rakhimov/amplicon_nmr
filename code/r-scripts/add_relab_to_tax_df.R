# The function add_relab_to_tax_df: 
# 1. Takes a table with absolute abundances
# 2. Groups rows by class and samples to calculate the number of reads in each 
# sample (TotalSample column is sum of Abundance within a sample). Group by class
# in case of repeat names in different classes
# 3. Then groups by class, sample, and a taxonomic rank (tax.rank) and
# calculates relative abundances: absolute abundances of each taxon in a sample 
# divided by the total number of reads in a sample, multiplied by 100. Save in a 
# RelativeAbundance column
# 4. Then, the dataframe is grouped by class and we create the TotalClass column:
# sum of reads from all samples belonging to a class.
# 5. Then, we group by class and the lowest taxonomic rank (tax.rank), and
# add a TotalAgglomRank column: sum of all reads of each tax.rank in each sample.
# 6. Finally, the MeanRelativeAbundance column is calculated as follows:
# MeanRelativeAbundance=TotalAgglomRank/TotalClass*100 
# sum of reads of a taxon in a host divided by sum of all reads from all taxa in a host 
# The sdRelativeAbundance is standard deviation of a taxon in a host
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
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100,
           sdRelativeAbundance=sd(RelativeAbundance))%>%
    ungroup()
  return(tax.df)
}
