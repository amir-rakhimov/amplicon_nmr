library(qiime2R)
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(pals)
library(patchwork)
library(Polychrome) 

truncationlvl<-"284-203"
# 313-229
# 284-229
# 284-203
# 273-203
# 273-198
# 265-203
# 265-198
trunclvls<-c("313-229",
             "284-229",
             "284-203",
             "273-203",
             "273-198",
             "265-203",
             "265-198")




for (truncationlvl in trunclvls){
  authorname<-"yasuda"
  directory<-paste0("./data/",authorname,"-data/")
  qiimedir<-paste0(directory,"qiime/")
  # agglom.rank<-"Family"
  agglom.rank<-"Genus"
  # Import qza files and convert them into a phyloseq object ####
  ps<-qza_to_phyloseq(
    # metadata = paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"),
    features = paste0(qiimedir,authorname,"-table-trimmed-dada2-",
                      truncationlvl,".qza"),
    taxonomy = paste0(qiimedir,authorname,"-taxonomy-trimmed-dada2-",
                      truncationlvl,".qza"),
    tree = paste0(qiimedir,authorname,"-rooted-tree-trimmed-dada2-",
                  truncationlvl,".qza")
  )
  # add custom metadata cause previous command loses a sample for some reason
  custom.md<-read.table(paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"),
                       header = T)
  colnames(custom.md)[1]<-"Sample"
  custom.md<-custom.md%>%column_to_rownames(var = "Sample")
  
  sample_data(ps)<-custom.md
  # View(ps@sam_data)
  # read_q2metadata(paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"))%>%View()
  
  # Create a dataframe of relative abundances ####
  ## Agglomerate by agglom.rank, transform into rel.ab, then ####
  # create a dataframe
  ps.agg.rel<-ps %>%
    #subset_samples(class=="NMR") %>% # choose only naked mole-rat
    tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
    subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
    transform_sample_counts(function(x)100* x / sum(x)) %>% # transform into percentages
    psmelt()  # transform the phyloseq object into an R dataframe
  head(ps.agg.rel) 
  
  
  
  ## rename d__Bacteria ####
  ps.agg.rel$Kingdom<-
    gsub("d__","",ps.agg.rel$Kingdom)
  
  # Total counts for each sample ####
  ps.total<-ps %>%
    tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
    subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
    psmelt() %>% # transform the phyloseq object into an R dataframe
    group_by(Sample) %>% # group by Sample id to perform operations on samples
    summarise(Abundance=sum(Abundance)) # sum abundances by sample id
  
  # Substitute NA values
  ps.agg.rel[is.na(ps.agg.rel)]<-"NonNA"
  
  ## Create a vector with taxonomic ranks ####
  phylum<-which(colnames(ps.agg.rel) =="Phylum") # identify the position of phylum column
  agglom.rank.col<-which(colnames(ps.agg.rel) ==agglom.rank) # identify the 
                                                # position of agglom.rank column
  taxcols<-phylum:agglom.rank.col # vector of taxonomic ranks
  
  # Replace empty taxa with Unclassified+previous taxon ####
  # we go through each taxonomic rank from phylum to agglom.rank
  # we substitute "uncultured" or "unclassified" and make it like 
  # Unclassified (Previous taxonomic rank)
  # foo<-ps.agg.rel
  # ps.agg.rel<-foo
  for (i in 1:nrow(ps.agg.rel)){
    for (j in taxcols){
      if(ps.agg.rel[i,j]=="NonNA"){
        ps.agg.rel[i,j]<-paste0("Unclassified"," (",ps.agg.rel[i,j-1],")") 
      }
      if (grepl("Unclassified",ps.agg.rel[i, j-1])){
        ps.agg.rel[i,j]<-ps.agg.rel[i,j-1]
      }
      if (ps.agg.rel[i,j]=="uncultured"){
        ps.agg.rel[i,j]<-paste0("Uncultured"," (",ps.agg.rel[i,j-1],")")
      }
      if (grepl("Uncultured",ps.agg.rel[i, j-1])){
        ps.agg.rel[i,j]<-ps.agg.rel[i,j-1]
      }
    }
  }
  
  # sanity check: is total relative abundance of each sample 100%?
  ps.agg.rel %>%
    group_by(Sample) %>% # Group by sample id
    mutate(Abundance = as.numeric(Abundance))%>% # transform Abundance column into numeric
    summarise(Abundance = sum(Abundance)) %>% # Sum all abundances
    pull(Abundance) %>% # get only a vector of numbers
    `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
    all() # get one TRUE/FALSE value
  
  # Extract taxa with highest mean relative abundance ####
  # dataset of three columns: two taxonomic ranks (e.g Genus, Family) 
  # and abundance 
  nmr.agg.rel<-ps.agg.rel %>%
    filter(class=="NMR")
  m.agg.rel<-ps.agg.rel %>%
    filter(class=="control")
  
  nmr.agg.rel.pooled<-nmr.agg.rel %>%
    group_by_at(c(agglom.rank.col,agglom.rank.col-1)) %>% # group by agglom.rank 
                                      # and the preceding column (based on index)
    mutate(Abundance = as.numeric(Abundance))%>% # convert Abundance column into numeric
    summarise(Abundance = mean(Abundance)) %>% # convert Abundance into mean abundance
    arrange(-Abundance) # order descending
  head(nmr.agg.rel.pooled)
  
  m.agg.rel.pooled<-m.agg.rel %>%
    group_by_at(c(agglom.rank.col,agglom.rank.col-1)) %>% # group by agglom.rank 
    # and the preceding column (based on index)
    mutate(Abundance = as.numeric(Abundance))%>% # convert Abundance column into numeric
    summarise(Abundance = mean(Abundance)) %>% # convert Abundance into mean abundance
    arrange(-Abundance) # order descending
  head(m.agg.rel.pooled)
  
  # get otus with mean abundance >1% ####
  nmr.1pc<-nmr.agg.rel.pooled[nmr.agg.rel.pooled$Abundance>1,]
  m.1pc<-m.agg.rel.pooled[m.agg.rel.pooled$Abundance>1,]
  
  # Create a Taxon column where items are either Unclassified (previous rank) ####
  # or agglom.rank (previous taxonomic rank)
  # e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
  nmr.agg.rel$Taxon<-ifelse(grepl("Unclassified|Uncultured",
                              nmr.agg.rel[[agglom.rank.col]]),  
                        nmr.agg.rel[[agglom.rank.col]], 
                        paste0(nmr.agg.rel[[agglom.rank.col]]," (",nmr.agg.rel[[agglom.rank.col-1]],")")) 
  # grepl finds which agglom.rank column items are uncultured or unclassified
  # if true, keeps Unclassified (previous rank)
  # if false (so, it's not Unclassified, actually has taxonomic classification), 
  # converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)
  m.agg.rel$Taxon<-ifelse(grepl("Unclassified|Uncultured",
                                 m.agg.rel[[agglom.rank.col]]),  
                           m.agg.rel[[agglom.rank.col]], 
                           paste0(m.agg.rel[[agglom.rank.col]]," (",m.agg.rel[[agglom.rank.col-1]],")")) 
  
  
  
  # If otu has mean rel.ab <1%, set it as Remainder ####
  nmr.agg.rel$Taxon<-ifelse(nmr.agg.rel[[agglom.rank.col]] %in% nmr.1pc[[1]],
                            nmr.agg.rel$Taxon,
                            "Remainder (Mean abundance < 1%)")
  m.agg.rel$Taxon<-ifelse(m.agg.rel[[agglom.rank.col]] %in% m.1pc[[1]],
                           m.agg.rel$Taxon,
                           "Remainder (Mean abundance < 1%)")
  # if our Taxon is in the 1pc dataset (mean abundance >1%), keep it as it is
  # otherwise, set it as Remainder (Mean abundance < 1%)
  
  
  # Create a dataframe of absolute abundances ####
  # remove taxa with 0 abundance
  ps.nonzero<-prune_taxa(taxa_sums(ps)>0,ps)
  
  ## extract absolute abundances ####
  ps.agg.abs<-ps.nonzero %>%
    #subset_samples(class=="NMR") %>% # choose only naked mole-rat
    tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
    subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
    psmelt()  # transform the phyloseq object into an R dataframe
  
  ps.agg.abs$Kingdom<-
    gsub("d__","",ps.agg.abs$Kingdom)
  
  ps.agg.abs[is.na(ps.agg.abs)]<-"NonNA"
  
  ## Replace empty taxa with Unclassified+previous taxon ####
  # find the column with your rank of interest
  # then extract columns between phylum and your rank
  phylum<-which(colnames(ps.agg.abs) =="Phylum")
  agglom.rank.col<-which(colnames(ps.agg.abs) ==agglom.rank)
  taxcols<-phylum:agglom.rank.col
  
  for (i in 1:nrow(ps.agg.abs)){
    for (j in taxcols){
      if(ps.agg.abs[i,j]=="NonNA"){
        ps.agg.abs[i,j]<-paste0("Unclassified"," (",ps.agg.abs[i,j-1],")") 
      }
      if (grepl("Unclassified",ps.agg.abs[i, j-1])){
        ps.agg.abs[i,j]<-ps.agg.abs[i,j-1]
      }
      if (ps.agg.abs[i,j]=="uncultured"){
        ps.agg.abs[i,j]<-paste0("Uncultured"," (",ps.agg.abs[i,j-1],")")
      }
      if (grepl("Uncultured",ps.agg.abs[i, j-1])){
        ps.agg.abs[i,j]<-ps.agg.abs[i,j-1]
      }
    }
  }
  
  ## add a column with merged taxonomic ranks ####
  ps.agg.abs$taxa.full<-
    apply(ps.agg.abs[,c((phylum-1):agglom.rank.col)],
          1,paste,collapse="_")
  
  save.image(paste0("./rdafiles/yasuda-",truncationlvl,"-",agglom.rank,
                    "-phyloseq-workspace.RData"))
  to.remove <- ls()
  to.remove <- c(to.remove[!grepl("^trunclvls", to.remove)], "to.remove")
  rm(list=to.remove)
}


