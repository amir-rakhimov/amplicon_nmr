# The function add_agegroup_to_tax_df takes a dataframe of abundances,
# the lowest taxonomic rank in the dataframe (e.g. genus), and a metadata with
# age groups as input. It joins the dataframe with metadata, groups by age 
# groups from the metadata, then sums the abundances from each age group.
# These sums are added as a column TotalAgegroup. Then, the function 
# groups by age group and the lowest taxonomic rank, and sums abundances for each
# taxonomic rank in each age group (TotalAgglomRankAge column). Finally,
# the function calculates the average relative abundance of each taxonomic rank
# inside each age group.
add_agegroup_to_tax_df<-function(tax.df,tax.rank,metadata.df){
  tax.df<-tax.df%>%
    left_join(metadata.df[,c("Sample","agegroup","old_agegroup")],by="Sample")%>%
    group_by(agegroup)%>% # group by class (animal host),
    mutate(TotalAgegroup=sum(Abundance))%>%
    group_by_at(c("agegroup",tax.rank))%>%
    mutate(TotalAgglomRankAge=sum(Abundance))%>%
    mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)%>%
    ungroup()%>%
    select(-TotalAgegroup,-TotalAgglomRankAge)
  return(tax.df)
}
