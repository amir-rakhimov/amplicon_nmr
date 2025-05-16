create_barplot_legend<-function(tax.df, # the dataset of abundances
                                inside.host, # whether you're comparing groups inside the host (Bool)
                                is.metagenome, # whether you're using metagenomic data (Bool)
                                split.by.comparison, # whether you want the legend to show unique taxa for compared groups
                                comparison.var.name, # what variable is used for grouping (class, age, sex)
                                tax.rank, # what taxonomic rank is used for the barplot
                                left.group.name, # name of the first group (inside comparison.var.name)
                                right.group.name, # name of the second group
                                taxa.to.plot,# list of taxa to be used in the legend
                                metadata.table # metadata for comparing groups inside host
                                ){ 
  # Left and right sets are for cases when want to highlight group-specific taxa on the barplot.
  # For example, when we compare NMR and non-NMR samples, or when we compare 
  # NMR age groups.
  # Find all unique taxa in left and right sets (even rare).
  if(inside.host==TRUE){
    # if we're comparing NMR groups (e.g. age)
    ### Find all unique taxa in left set samples ####
    left.set<-tax.df%>%
      left_join(metadata.table)%>%
      ungroup()%>%
      filter(get(comparison.var.name)==left.group.name)%>%
      select(all_of(tax.rank))%>%
      distinct()%>%
      pull
    ### Find all unique taxa in right set samples ####
    right.set<-tax.df%>%
      left_join(metadata.table)%>%
      ungroup()%>%
      filter(get(comparison.var.name)==right.group.name)%>%
      select(all_of(tax.rank))%>%
      distinct()%>%
      pull
  }else{
    ### Find all unique taxa in left set samples (not inside host) ####
    left.set<-tax.df%>%
      filter(class==left.group.name)%>%
      select(all_of(tax.rank))%>%
      unique()%>%
      pull()
    ### Find all unique taxa in right set samples (not inside host) ####
    right.set<-tax.df%>%
      filter(class!=left.group.name)%>%
      select(all_of(tax.rank))%>%
      unique()%>%
      pull()
  }
  ### 8.3  Get NMR-specific taxa ####
  left.set.uniq<-setdiff(left.set,right.set)
  right.set.uniq<-setdiff(right.set,left.set)
  
  #  Obtain taxa found in the taxa.to.plot vector (sorted unclassified and  ####
  # classified taxa) with MeanRelativeAbundance>=1%. Keep only unique values in the
  # tax.rank.vec vector
  tax.rank.vec<-tax.df%>%
    filter(Taxon.bp%in%taxa.to.plot,MeanRelativeAbundance>=1)%>%
    select(tax.rank)%>%
    pull()%>%
    unique()
  
  ### 9.1 Find NMR-specific taxa in the tax.rank.vec vector  ####
  left.set.uniq.legend<-tax.rank.vec[tax.rank.vec%in%left.set.uniq]
  left.set.uniq.legend<-tax.df%>%
    filter(Taxon.bp%in%taxa.to.plot,MeanRelativeAbundance>=1)%>%
    distinct(get(tax.rank),Taxon.bp)%>%
    rename(!!tax.rank:="get(tax.rank)")%>%
    filter(get(tax.rank)%in%left.set.uniq)%>%
    select(Taxon.bp)%>%
    pull()
  right.set.uniq.legend<-tax.rank.vec[tax.rank.vec%in%right.set.uniq]
  right.set.uniq.legend<-tax.df%>%
    filter(Taxon.bp%in%taxa.to.plot,MeanRelativeAbundance>=1)%>%
    distinct(get(tax.rank),Taxon.bp)%>%
    rename(!!tax.rank:="get(tax.rank)")%>%
    filter(get(tax.rank)%in%right.set.uniq)%>%
    select(Taxon.bp)%>%
    pull()
  
  final.legend<-as.data.frame(taxa.to.plot)%>%
    mutate(new.colors=taxa.to.plot)
  if(is.metagenome==TRUE){
    # If data is metagenomic
    final.legend<-final.legend%>%
      left_join(unique(subset(tax.df[,c("Kingdom","Family","Genus","Species","Taxon.bp")],
                              Taxon.bp%in%taxa.to.plot[-1])),
                by=join_by(new.colors==Taxon.bp))%>%
      # mutate(new.colors=Taxon.bp)%>%
      mutate(Taxon.bp=new.colors)%>%
      mutate(new.colors=ifelse(Kingdom=="Bacteria",
                               paste0("<span style='color: orange'><b><i>",Taxon.bp,"</i></b></span>"),
                               new.colors),
             new.colors=ifelse(Kingdom=="Eukaryota",
                               paste0("<span style='color: palegreen4'><b><i>",Taxon.bp,"</i></b></span>"),
                               new.colors),
             new.colors=ifelse(Kingdom=="Archaea",
                               paste0("<span style='color: violetred4'><b><i>",Taxon.bp,"</i></b></span>"),
                               new.colors))%>%
      select(Kingdom,Species,Genus,Taxon.bp,new.colors)
    final.legend[1,]<-data.frame(Kingdom="Remainder (Mean relative abundance < 1%)",
                                 Species="Remainder (Mean relative abundance < 1%)",
                                 Genus="Remainder (Mean relative abundance < 1%)",
                                 Taxon.bp="Remainder (Mean relative abundance < 1%)",
                                 new.colors="Remainder (Mean relative abundance < 1%)")
    
  }else{
    # If data isn't metagenomic
    final.legend<-final.legend%>%
      rename("Taxon.bp"="taxa.to.plot")
  }
  if(split.by.comparison==TRUE){
    # If you want to show unique taxa for compared groups
    final.legend<-final.legend%>%
      mutate(new.colors=ifelse(Taxon.bp%in%left.set.uniq.legend,
                               paste0("<span style='color: firebrick'><b><i>",Taxon.bp,"</i></b></span>"),
                               new.colors),
             new.colors=ifelse(Taxon.bp%in%right.set.uniq.legend&right.group.name!="",
                               paste0("<span style='color: dodgerblue4'><b><i>",Taxon.bp,"</i></b></span>"),
                               new.colors))
  }
  final.legend<-final.legend%>%
    mutate(new.colors=gsub("_"," ", new.colors))
  
  return(final.legend)
}