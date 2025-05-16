# The function create_summary_stats_table takes a dataframe
# with ASV abundances (tax.df), groups by class (animal host), and adds a column with 
# number of samples per host (TotalSamplesPerHost). It uses a dplyr function n_distinct
# 2. Then, the function adds a column with number of total reads per host (TotalReadsPerHost).
# It uses a function sum() on Abundance column
# 3. Then, the function groups by Sample and adds a column with number of total reads 
# per sample (LibrarySize). It uses sum() on Abundance column
# 4. Then, the function keeps unique samples with distinct() function
# 5. Then, the function groups by class and adds a column with average number of reads
# per sample (MeanLibrarySize). It uses mean() function on LibrarySize.
# 6. The function also adds a Standard deviation of library size (SDLibrarySize). The 
# rationale is the same as when calculating MeanLibrarySize, except the function is
# sd()
# 7. The function selects columns class, TotalSamplesPerHost, TotalReadsPerHost, 
# MeanLibrarySize, and SDLibrarySize.
# 8. The function selects unique rows (distinct() function).
# 9. The function sorts by class column, then adds several dataframes:
# 9.1: Dataframe with the number of unique ASVs per host (n.asv.table).
# It is added as `n` column, but is renamed to ASVPerHost.
# 9.2: Dataframe with the number of unique phyla per host (n.phylum.table).
# It is added as `n` column, but is renamed to PhylaPerHost.
# 9.3: Dataframe with the number of unique families per host (n.family.table).
# It is added as `n` column, but is renamed to FamiliesPerHost.
# 9.4: Dataframe with the number of unique genera per host (n.genus.table).
# It is added as `n` column, but is renamed to GeneraPerHost.
create_summary_stats_table<-function(tax.df,
                                     n.asv.table,
                                     n.phylum.table,
                                     n.family.table,
                                     n.genus.table){
  final.summary.stats.table<-tax.df%>%
    group_by(class)%>%
    mutate(TotalSamplesPerHost=n_distinct(Sample))%>%
    mutate(TotalReadsPerHost=sum(Abundance))%>%
    group_by(Sample)%>%
    mutate(LibrarySize=sum(Abundance))%>%
    distinct(Sample,.keep_all = T)%>%
    group_by(class)%>%
    mutate(MeanLibrarySize =round(mean(LibrarySize)),
           SDLibrarySize=round(sd(LibrarySize)))%>%
    select(class,
           TotalSamplesPerHost,
           TotalReadsPerHost,
           MeanLibrarySize,
           SDLibrarySize)%>%
    distinct(class,.keep_all = T)%>%
    arrange(class)%>%
    left_join(n.asv.table)%>%
    rename(ASVPerHost=n)%>%
    left_join(n.phylum.table)%>%
    rename(PhylaPerHost=n)%>%
    left_join(n.family.table)%>%
    rename(FamiliesPerHost=n)%>%
    left_join(n.genus.table)%>%
    rename(GeneraPerHost=n)
  return(final.summary.stats.table)
}
