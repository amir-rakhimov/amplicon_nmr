# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan"))
library(tidyverse)
library(phyloseq)
library(vegan)

truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
authorname<-"pooled"
date_time<-"20240426_21_44_30"

# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
load(file.path("./output/rdafiles",paste(
  date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

# custom.levels<-c("NMR","B6mouse")
pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "B6mouse" = "B6 mouse",
                       "MSMmouse" = "MSM/Ms mouse",
                       "FVBNmouse" = "FVB/N <br>mouse",
                       "DMR" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus <br>cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pvo" = "*Pteromys volans orii*",
                       "NMRwt"="Wild *Heterocephalus glaber*"
)
custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)

# If we're working with NMR only
# custom.levels<-c("NMR")
# custom.levels<-c("B6mouse",
#                  "MSMmouse",
#                  "FVBNmouse")

# Import data ####
ps.q.df <-ps.q.agg%>%
  select(Sample,Abundance,class,all_of(agglom.rank))%>%
  filter(Abundance!=0)

ps.q.df<-ps.q.df%>%
  filter(class %in% custom.levels,Abundance!=0)


# Convert the data frame into wide format: rows are samples and columns
# are taxa
ps.q.df.wide<-ps.q.df%>%
  pivot_wider(names_from = all_of(agglom.rank), # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# Set sample names as row names
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  


# Rarefaction ####
# find the smallest sample size
min.n_seqs.all<-ps.q.agg%>%
  filter(class %in% custom.levels)%>%
  select(Sample,all_of(agglom.rank),Abundance)%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# If you want to plot a rarefaction curve ####
# library(ggrepel)
# ps.q.mat<-as(t(otu_table(ps.q)),"matrix") # from phyloseq
# ps.q.mat<-as.matrix(ps.q.df.wide) # from the ps.q.df matrix
# set.seed(1)
# rare.df<-rarecurve(ps.q.mat,step = 100,sample=min(rowSums(ps.q.mat)),tidy = TRUE)
# rare.df%>%
#   filter(Sample<=100000)%>%
#   group_by(Site)%>%
#   mutate(label=if_else(Sample==max(Sample),as.character(Site),NA_character_))%>%
#   filter(Site%in%rownames(custom.md[which(custom.md$class=="NMR"),]))%>%
#   ggplot(.,aes(x=Sample,y=Species,col=Site))+
#   geom_line()+
#   # coord_cartesian(xlim=c(0,100000))+
#   geom_vline(xintercept = min(rowSums(ps.q.mat)))+
#   annotate("text",
#            x=min(rowSums(ps.q.mat))+2000,
#            y=10,
#            label=min(rowSums(ps.q.mat)))+
#   geom_label_repel(aes(label = label),
#                    nudge_x = 1,
#                    na.rm = TRUE) +
#   theme_bw()+
#   labs(x="Sample size",
#        y="ASV")+
#   theme(legend.position = "none")
# ggsave(paste0("./images/lineplots/",
#               paste(paste(format(Sys.time(),format="%Y%m%d"),
#                           format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                     "rarecurve",
#                     truncationlvl,agglom.rank,
#                     sep = "-"),".png"),
#        plot=last_plot(),
#        width = 4500,height = 3000,
#        units = "px",dpi=300,device = "png")

# rarefied asv table with vegan
set.seed(1)
ps.q.df.rare<-rrarefy(ps.q.df.wide,sample=min.n_seqs.all)
ps.q.df.rare<-ps.q.df.rare%>%
  as_tibble(rownames="Sample")%>%
  pivot_longer(-Sample)%>%
  as.data.frame()%>%
  left_join(unique(ps.q.agg[,c("Sample","class","sex","birthday")]),
             by="Sample")
if(asvlevel==TRUE){
  ps.q.df.rare<-ps.q.df.rare%>%
    rename(OTU=name,Abundance=value)%>%
    filter(Abundance!=0)  
}else{
  # rename the 'name' column corresponding to the agglom.rank
  ps.q.df.rare[,paste(agglom.rank)]<-ps.q.df.rare$name
  ps.q.df.rare<-ps.q.df.rare%>%
    select(-name)%>%
    rename(Abundance=value)%>%
    filter(Abundance!=0)
}
write.table(ps.q.df.rare,
            file = file.path("./output/rtables",authorname,paste0(
              paste(
                paste(format(Sys.time(),format="%Y%m%d"),
                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                "ps.q.df.rare-nonfiltered",agglom.rank,
                paste(custom.levels,collapse = '-'),sep = "-"),
              ".tsv")),
            row.names = F,
            sep = "\t")


# Filtering by prevalence ####
# for each host, calculate the number of samples an ASV was observed in 
# so the same ASV can be observed in multiple hosts.
# observed_samples shows number of samples per host
ps.q.sample_counts <- ps.q.agg%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class,Family,Genus) %>% 
  summarize(observed_samples = n_distinct(Sample))

# calculate the total number of samples for each host
ps.q.host_counts <- ps.q.agg%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class)%>%
  summarize(total_samples = n_distinct(Sample)) %>%
  ungroup()

# Join host counts with sample counts
# Create the PercentageSamples column: num of samples  from a specific host
# an ASV was observed in divided by the total number of samples in that host
# Then, we join our dataframe with mean relative abundance df (ps.agg.rel.pooled).
# Then, we add the Taxon column from ps.q.agg.abs
ps.q.prevalences <- ps.q.host_counts %>%
  left_join(ps.q.sample_counts, by = "class") %>%
  mutate(PercentageSamples = observed_samples/total_samples*100)%>% # percentage
  # of samples an ASV was observed in
  left_join(ps.q.agg[,c("class","Genus","Family","Taxon","MeanRelativeAbundance")],
            by=c("class","Genus","Family")) # add mean
  # relative abundances and the Taxon column

## filter by percentage of samples an ASV was observed in ####
ps.q.df.norare.filtered <-ps.q.df%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"),relationship = "many-to-many")%>%
  filter(PercentageSamples>=10)
ps.q.df.rare.filtered<-ps.q.df.rare%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"),relationship = "many-to-many")%>%
  filter(PercentageSamples>=10,MeanRelativeAbundance>=1)


write.table(ps.q.df.norare.filtered,
            file = paste0("./rtables/",authorname,"/ps.q.df.norare.filtered-",
                          paste(agglom.rank,custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            sep = "\t")
write.table(ps.q.df.rare.filtered,
            file = paste0("./rtables/",authorname,"/ps.q.df.rare.filtered-",
                          paste(agglom.rank,custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            sep = "\t")
