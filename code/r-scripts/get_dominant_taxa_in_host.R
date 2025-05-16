# The function get_dominant_taxa_in_host uses a dataframe with abundances, 
# the taxonomic rank for which you want to find dominant taxa (e.g. dominant genera), 
# and the host name/names.
# If there is on host, the function filters the dataframe, adds a column
# with the total abundance (sum of all reads in the host), then groups
# by taxonomic rank (e.g. genus) and calculates the total abundance of each taxon.
# After that, the function adds a column with average relative abundance 
# of each taxon.
# If there are two or more hosts, the procedure is the same, except we 
# group by class before counting total abundance and group by both class and 
# taxonomic rank before calculating total abundance of each taxon.
# The function selects the necessary columns like Phylum or Family, sorts
# by average relative abundance and returns a filtered dataframe.
get_dominant_taxa_in_host<-function(tax.df,tax.rank,host,nonbacterial.table){
  if(length(host)==1){
    dominant_taxa_table<-tax.df%>%
      ungroup()%>%
      filter(class==host)%>%
      mutate(TotalClass=sum(Abundance))%>%
      group_by_at(c(tax.rank))%>%
      mutate(TotalAgglomRank=sum(Abundance))%>%
      mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
    
  }else if (length(host)>1){
    dominant_taxa_table<-tax.df%>%
      filter(class%in%host)%>%
      group_by(class)%>%
      mutate(TotalClass=sum(Abundance))%>%
      group_by_at(c("class",tax.rank))%>%
      mutate(TotalAgglomRank=sum(Abundance))%>%
      mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  }
  if(nonbacterial.table==TRUE){
    dominant_taxa_table<-dominant_taxa_table%>%
      filter(Kingdom!="Bacteria")
  }
  
  if(tax.rank=="Phylum"){
    dominant_taxa_table<-dominant_taxa_table%>%
      select(Phylum,MeanRelativeAbundance)
  }else if (tax.rank=="Family"){
    dominant_taxa_table<-dominant_taxa_table%>%
      select(Phylum,Family,MeanRelativeAbundance)
  }else if (tax.rank=="Genus"){
    dominant_taxa_table<-dominant_taxa_table%>%
      select(Phylum,Family,Genus,MeanRelativeAbundance)
  }else if (tax.rank=="Species"){
    dominant_taxa_table<-dominant_taxa_table%>%
      select(Species,MeanRelativeAbundance)
  }
  dominant_taxa_table<-dominant_taxa_table%>%
    distinct()%>%
    arrange(-MeanRelativeAbundance)
  return(dominant_taxa_table)
}
