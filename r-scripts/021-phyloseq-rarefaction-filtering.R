library(tidyverse)
library(phyloseq)
library(vegan)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

custom.levels<-c("NMR","SPFmouse")
custom.levels<-c("NMR","SPFmouse","spalax","FukomysDamarensis","rabbit")

# Import data ####
ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)%>%
  filter(Abundance!=0)

ps.q.df<-ps.q.df%>%
  filter(class %in% custom.levels)


# convert the data frame into wide format
ps.q.df.wide<-ps.q.df%>%
  select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  


# Rarefaction ####
# find the smallest sample size
min.n_seqs.all<-ps.q.agg.abs%>%
  filter(class %in% custom.levels)%>%
  select(Sample,OTU,Abundance)%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# rarefied asv table with vegan
set.seed(1)
ps.q.df.rare<-rrarefy(ps.q.df.wide,sample=min.n_seqs.all)%>%
  as_tibble(rownames="Sample")%>%
  pivot_longer(-Sample)%>%
  as.data.frame()%>%
  inner_join(unique(ps.q.agg.abs[,c("Sample","class")]),
             by="Sample")%>%
  rename(Taxon=name,Abundance=value)

# Filtering by prevalence ####
# for each host, calculate the number of samples an ASV was observed in 
ps.q.sample_counts <- ps.q.agg.abs%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class,Family,Genus) %>% 
  summarize(observed_samples = n_distinct(Sample))

# calculate the total number of samples for each host
ps.q.host_counts <- ps.q.agg.abs%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class)%>%
  summarize(total_samples = n_distinct(Sample)) %>%
  ungroup()

ps.q.prevalences <- ps.q.host_counts %>%
  left_join(ps.q.sample_counts, by = "class") %>%
  mutate(PercentageSamples = observed_samples/total_samples*100)%>% # percentage
  # of samples an ASV was observed in
  left_join(ps.q.agg.rel.pooled,by=c("class","Genus","Family"))%>% # add mean
  # relative abundances
  rename(MeanAbundance=Abundance)%>%
  left_join(unique(ps.q.agg.abs[,c("Genus","Family","Taxon","class")]),
            by=c("class","Genus","Family")) # add the Taxon column

## filter by percentage of samples an ASV was observed in ####
ps.q.df.norare.filtered <-ps.q.df%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"))%>%
  filter(PercentageSamples>=10)
ps.q.df.rare.filtered<-ps.q.df.rare%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"))%>%
  filter(PercentageSamples>=10)


write.table(ps.q.df.norare.filtered,
            file = paste0("./rtables/alldir/ps.q.df.norare.filtered-",
                          paste(custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            sep = "\t")
write.table(ps.q.df.rare.filtered,
            file = paste0("./rtables/alldir/ps.q.df.rare.filtered-",
                          paste(custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            sep = "\t")
