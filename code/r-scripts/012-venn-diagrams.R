library(phyloseq)
library(tidyverse)
asvlevel=FALSE
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "B6mouse" = "B6 mouse",
                       "MSMmouse" = "MSM/Ms mouse",
                       "FVBNmouse" = "FVB/N mouse",
                       "DMR" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pvo" = "*Pteromys volans orii*")
custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)

ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)%>%
  filter(Abundance!=0)

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
  summarize(total_samples = n_distinct(Sample)) %>% # num of total_samples for each host
  ungroup()

ps.q.prevalences <- ps.q.host_counts %>% # take the num of samples for each host
  left_join(ps.q.sample_counts, by = "class") %>% # join with taxa and 
                              # num of samples each taxon was observed in
                            # each host
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
  filter(PercentageSamples>=10)%>%
  filter(MeanAbundance>=1)

taxa.list<-list()
for(i in seq_along(custom.levels)){ # for each host
  lvl.df<-ps.q.agg.rel%>% # take the dataset of relative abundances
    left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data 
                                         # with the dataset of total abundances, 
                                        # so we can have info about sample size
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change 
                                            # the order of our class column, 
                                                # so the NMR will be first
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of 
    # our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")"))%>% # add a
                                        # column with sample sizes for each ASV
    filter(class==custom.levels[i],Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[i])
  lvl.name<-gsub("<br>"," ", lvl.name)
  taxa.list[[custom.levels[i]]]<-unique(lvl.df$Taxon) # get unique taxa
}

for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.df.norare.filtered%>% # take the filtered dataset
    filter(class==custom.levels[i],Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[i])
  lvl.name<-gsub("<br>"," ", lvl.name)
  taxa.list[[custom.levels[i]]]<-unique(lvl.df$Taxon)
}

# vennD<-ggVennDiagram(taxa.list, 
#                      label_alpha = 0,
#                      category.names = pretty.facet.labels, 
#                      label_geom = "label", 
#                      label_size = 4) +
#   scale_fill_distiller(palette = "Reds", direction = 1)+ 
#   scale_x_continuous(expand = expansion(mult = .2))
# 
# ggsave(paste0("./images/venn/",Sys.Date(),
#               "-venn-",paste(names(taxa.list),collapse = "_"),"-",
#               agglom.rank,"-",truncationlvl,".png"),plot = vennD,
#        width = 2000,height = 2500,
#        dpi = 300,units = "px",device = "png") 

# Find unique taxa
foo<-ps.q.agg.abs
ps.q.agg.abs<-foo

ps.q.agg.abs<-ps.q.agg.abs%>%
  filter(class%in%c("B6mouse","MSMmouse","FVBNmouse"))

ref.class<-"MSMmouse" # find unique taxa for this class
ref.class<-"B6mouse" # find unique taxa for this class
ref.class<-"FVBNmouse" # find unique taxa for this class
ref.class<-"NMR" # find unique taxa for this class

taxa.no_ref<-ps.q.agg.abs%>%
  filter(class!=ref.class,Abundance!=0)%>%
  select(Taxon)%>%
  distinct()%>%
  as_vector()
taxa.ref<-ps.q.agg.abs%>%
  filter(class==ref.class,Abundance!=0)%>%
  select(Taxon)%>%
  distinct()%>%
  as_vector()
uniq.taxa<-setdiff(taxa.ref,taxa.no_ref) # Unique for ref.class
intersect(setdiff(taxa.ref,taxa.no_ref),ps.q.df.norare.filtered$Taxon)
setdiff(taxa.no_ref,taxa.ref)

ps.q.agg.abs[ps.q.agg.abs$Taxon%in%uniq.taxa,]%>%
  filter(class==ref.class,Abundance!=0)%>%View

# Find unique taxa in filtered data
taxa.no_nmr<-ps.q.df.norare.filtered%>%
  filter(class!=ref.class,Abundance!=0)%>%
  select(Taxon)%>%
  distinct()%>%
  as_vector()
taxa.nmr<-ps.q.df.norare.filtered%>%
  filter(class==ref.class,Abundance!=0)%>%
  select(Taxon)%>%
  distinct()%>%
  as_vector()
foo<-setdiff(taxa.nmr,taxa.no_nmr)
ps.q.agg.abs[ps.q.agg.abs$Taxon%in%foo,]%>%
  filter(Abundance!=0)%>%View


# Find unique asv
if (asvlevel==TRUE){
  asv.no_nmr<-ps.q.agg.abs%>%
    filter(class!=ref.class,Abundance!=0)%>%
    select(OTU)%>%
    distinct()%>%
    as_vector()
  asv.nmr<-ps.q.agg.abs%>%
    filter(class==ref.class,Abundance!=0)%>%
    select(OTU)%>%
    distinct()%>%
    as_vector()
  setdiff(asv.nmr,asv.no_nmr) # No unique asv
  setdiff(asv.no_nmr,asv.nmr) 
}


