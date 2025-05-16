# In filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank], collapse='|'),get(agglom.rank))),
# we are searching agglom.rank column for taxa that have a higher rank
# in their name. For example, Bacteria Kingdom is unclassified, and if it's in
# the agglom.rank column, we find the word "Bacteria" which makes it an 
# unclassified taxon.
# In summarise(TotalUnclassifiedPercent=sum(RelativeAbundance)), we sum the 
# relative abundance of all taxa that were unclassified in a given sample
# from a given host (because we group by sample and class).
# We calculate Mean, SD, min, max, and median of unclassified percentages.
# After that, there is no need for the TotalUnclassifiedPercent column because 
# the summary statistics are calculated for every host. So, we keep unique rows
# Note: stats are the same for rarefied and non-rarefied data
get_unclassified_summary_stats<-function(tax.df,
                                         tax.rank){
  all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
  # % of unclassified data in each sample
  unclassified.percent.by.sample<-tax.df%>%
    filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                       collapse='|'),get(agglom.rank)))%>%
    group_by(Sample,class)%>%
    summarise(TotalUnclassifiedPercent=sum(RelativeAbundance))%>%
    ungroup()
  # Stats of unclassified % (min, max, mean, SD, median)
  unclassified.stats.by.class<-unclassified.percent.by.sample%>%
    group_by(class)%>%
    summarise(MeanTotalUnclassifiedPercent=round(mean(TotalUnclassifiedPercent)),
           SDTotalUnclassifiedPercent=round(sd(TotalUnclassifiedPercent)),
           minTotalUnclassifiedPercent=round(min(TotalUnclassifiedPercent)),
           maxTotalUnclassifiedPercent=round(max(TotalUnclassifiedPercent)),
           MedianTotalUnclassifiedPercent=round(median(TotalUnclassifiedPercent)))%>%
    # select(-Sample,-TotalUnclassifiedPercent)%>%
    # distinct(class,.keep_all = T)%>%
    arrange(-MeanTotalUnclassifiedPercent)
  
  ### Get the number of unclassified taxa in each host ####
  unclassified.taxa.n<-tax.df%>%
    filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                       collapse='|'),get(agglom.rank)))%>%
    group_by(class)%>%
    distinct(get(tax.rank),.keep_all = T)%>%
    tally()%>%
    arrange(-n)%>%
    rename(NumUnclassifiedTaxa=n)
  
  ### Get the number of classified taxa in each host (unrarefied) ####
  classified.taxa.n<-tax.df%>%
    filter(!grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                        collapse='|'),get(agglom.rank)))%>%
    group_by(class)%>%
    distinct(get(tax.rank),.keep_all = T)%>%
    tally()%>%
    arrange(-n)%>%
    rename(NumClassifiedTaxa=n)
  
  total.taxa.by.class<-tax.df%>%
    group_by(Sample,class)%>%
    distinct(get(tax.rank),.keep_all = T)%>%
    tally()%>%
    group_by(class)%>%
    summarise(MinTotalTaxa=min(n),
              MaxTotalTaxa=max(n))
  
  unclassified.taxa.summary.stats.table<-unclassified.stats.by.class%>%
    left_join(unclassified.taxa.n)%>%
    left_join(classified.taxa.n)%>%
    left_join(total.taxa.by.class)
  
  
  return(unclassified.taxa.summary.stats.table)
}
