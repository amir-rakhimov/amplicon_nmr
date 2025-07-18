# Compare our data with results from the original papers
library(tidyverse)
library(phyloseq)
# library(DT)
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/create_summary_stats_table.R")
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
source("./code/r-scripts/add_relab_to_tax_df.R")
source("./code/r-scripts/add_agegroup_to_tax_df.R")
source("./code/r-scripts/get_unclassified_summary_stats.R")
source("./code/r-scripts/ggplot_species.R")
source("./code/r-scripts/add_zero_rows.R")
# Import my dataset ####
## 2. Specifying parameters and directory/file names #### 
rdafiles.directory<-"./output/rdafiles"

truncationlvl<-"234" # truncation level that we chose in QIIME2
# truncationlvl<-"0" # truncation level that we chose in QIIME2

authorname<-"pooled" # name of the folder with QIIME2 output
# authorname<-"biagi" # name of the folder with QIIME2 output

rtables.directory<-file.path("./output/rtables",authorname)
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
# Import datasets as rds files
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))


# Re-analyse the DMR dataset ####
dmr.asvs <- read_tsv("./FreezeDriedVSFrozen-main/data/asv_table_FreezedriedVsFrozen.tsv", 
                     col_types =cols(asv = col_character(),
                                     sample = col_character(),
                                     count = col_double(),
                                     relab = col_double() ))

dmr.metadata <- read_csv("./data/dmr-original-data/FDvsFrozenMetadata.csv", 
                         col_types = cols(.default = col_character(),
                                          SampleDate = col_date(format = ""),
                                          SampleOrder = col_number())) %>% 
  suppressWarnings()

# Rename Sample number
dmr.metadata  <- dmr.metadata %>% 
  dplyr::rename(OldSampleNumber = SampleNumber, 
                SampleNumber = NewSampleNumber )

# make a character of the sample number
dmr.metadata <- dmr.metadata %>% 
  mutate(Asample = paste("A", sample, sep= "_"))

dmr.taxonomy <- read_tsv("./data/dmr-original-data/taxonomy_FreezedriedVsFrozen.tsv", 
                         col_types =  cols(asv = col_character(),
                                           kingdom = col_character(), 
                                           phylum = col_character(),
                                           class = col_character(),
                                           order = col_character(),
                                           family = col_character(),
                                           genus = col_character(),
                                           species = col_character()))

# Rename unclassified
dmr.taxonomy  <- dmr.taxonomy %>%
  mutate(
    phylum = ifelse(is.na(phylum), sprintf("%s unclassified", kingdom), phylum),
    class = ifelse(is.na(class), sprintf("%s unclassified", str_remove(phylum, "unclassified")), class),
    order = ifelse(is.na(order), sprintf("%s unclassified", str_remove(class, "unclassified")), order),
    family = ifelse(is.na(family), sprintf("%s unclassified", str_remove(order, "unclassified")), family),
    genus = ifelse(is.na(genus), sprintf("%s unclassified", str_remove(family, "unclassified")), genus),
    species = ifelse(is.na(species), sprintf("%s unclassified", str_remove(genus, "unclassified")), species))





## Filter to keep only frozen samples ####
# Also, join the asv data with taxonomy table
dmr.asv.tax.df<-dmr.asvs%>%
  left_join(dmr.taxonomy)%>%
  rename("Sample"=sample)%>%
  rename("Abundance"=count)%>%
  rename("Kingdom"=kingdom)%>%
  rename("Phylum"=phylum)%>%
  rename("Class"=class)%>%
  rename("Order"=order)%>%
  rename("Family"=family)%>%
  rename("Genus"=genus)%>%
  rename("Species"=species)%>%
  mutate(class="DMR")

dmr.asv.tax.df.frozen<-dmr.asv.tax.df%>%
  filter(Sample%in%subset(dmr.metadata,Treatment=="Frozen")$sample)

dmr.taxonomy.frozen<-dmr.taxonomy%>%
  filter(asv%in%dmr.asv.tax.df.frozen$asv)

## ====== Results after Table 1: ====== ####
## Compare numbers of phyla, families, genera, and asv between datasets ####
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Phylum")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Phylum")

get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Family")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Family")

get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Genus")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Genus")

get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "asv")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "asv")

## Classified and unclassified genera ####
# Unclassified
dmr.taxonomy.frozen[grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 41 unclassified genus in frozen

dmr.taxonomy[grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 51 unclassified genus in all


# Classified
dmr.taxonomy.frozen[!grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 124 classified genus in frozen
dmr.taxonomy[!grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 159 classified genus in all


## How many total reads in frozen data ####
# 1 853 295 reads in frozen data
sum(dmr.asv.tax.df.frozen$Abundance)
# 3 626 584 reads in total
sum(dmr.asv.tax.df$Abundance)

## Library size ####
# Mean 92665, SD 44202
dmr.asv.tax.df.frozen%>%
  group_by(Sample)%>%
  summarise(LibrarySize=sum(Abundance))%>%
  summarise(MeanLibrarySize=mean(LibrarySize),
            SDLibrarySize=sd(LibrarySize))


## Dominant phyla in the frozen dataset (Results). Compare with my data ####
dmr.asv.tax.df.frozen %>% 
  group_by(Sample, Phylum) %>% 
  summarise(sample_relab = sum(relab)*100) %>%
  group_by(Phylum) %>% 
  summarise(mean_relab =mean(sample_relab), 
            min_relab= min(sample_relab), 
            max_relab= max(sample_relab))%>%
  arrange(desc(mean_relab)) %>% 
  mutate(mean_relab = round(mean_relab,digits = 5),
         min_relab = round(min_relab,digits = 5),
         max_relab = round(max_relab,digits = 5))

### Dominant phyla in my data (re-analysis in the Results) ####
ps.q.agg.phylum%>%
  filter(class=="DMR")%>% 
  group_by(Phylum) %>% 
  summarise(mean_relab =mean(RelativeAbundance ), 
            min_relab= min(RelativeAbundance ), 
            max_relab= max(RelativeAbundance ))%>%
  arrange(desc(mean_relab)) %>% 
  mutate(mean_relab = round(mean_relab,digits = 5),
         min_relab = round(min_relab,digits = 5),
         max_relab = round(max_relab,digits = 5))




# Re-analyse the PVO dataset ####
## Dominant phyla (Results section) ####
get_dominant_taxa_in_host(ps.q.agg.phylum,
                          tax.rank = "Phylum",
                          host = "pvo",
                          nonbacterial.table = F)
ps.q.agg.phylum%>%
  filter(class=="pvo")%>%
  ungroup()%>%
  group_by(Phylum)%>%
  summarise(mean_relab =mean(RelativeAbundance ), 
            min_relab= min(RelativeAbundance ), 
            max_relab= max(RelativeAbundance ))%>%
  arrange(desc(mean_relab)) %>% 
  mutate(mean_relab = round(mean_relab,digits = 5),
         min_relab = round(min_relab,digits = 5),
         max_relab = round(max_relab,digits = 5))

## Dominant families ####
get_dominant_taxa_in_host(ps.q.agg.family,
                          tax.rank = "Family",
                          host = "pvo",
                          nonbacterial.table = F)


# Re-analyse the Hare and rabbits dataset ####
custom.md%>%
  filter(class%in%c("hare","rabbit"))%>%
  arrange(class,Sample)
## Remove pregnant and lactating samples
selected.hare_rabbit.samples<-c("Hare 1","Hare 2", "Hare 3","Hare 4","Hare 5",
                                "Hare 6","Hare 7","Hare 8",
                                "Rabbit 2","Rabbit 3", "Rabbit 6","Rabbit 7",
                                "Rabbit 10", "Rabbit 11", "Rabbit 12")

## Import abundances from the original paper (supplementary data). Only 
# Illumina samples
hare_rabbit.phylum<-read.table("./data/hares-rabbits-original-data/phyla.tsv",sep = "\t",
                               header = T)%>%
  as_tibble()
# Change to long format and remove the "Unassigned" taxon
hare_rabbit.phylum<-hare_rabbit.phylum%>%
  filter(Sample%in%selected.hare_rabbit.samples)%>%
  pivot_longer(!Sample,names_to = "Phylum",values_to = "Abundance")%>%
  filter(Phylum!="Unassigned")%>%
  filter(Abundance!=0)

### Import family abundances ####
hare_rabbit.family<-read.table("./data/hares-rabbits-original-data/families.tsv",
                               sep = "\t",
                               header = T)%>%
  as_tibble()
# Remove rows with zeros
hare_rabbit.family<-hare_rabbit.family%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Family",values_to = "Abundance")%>%
  filter(Family!="unassigned")%>%
  filter(Abundance!=0)

### Import genus abundances ####
hare_rabbit.genus<-read.table("./data/hares-rabbits-original-data/genera.tsv",sep = "\t",
                              header = T)%>%
  as_tibble()

# Remove rows with zeros
hare_rabbit.genus<-hare_rabbit.genus%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Genus",values_to = "Abundance")%>%
  filter(Abundance!=0)
 

## Count taxa in the original data vs my re-analysis ####
# 14 phyla in the original data
hare_rabbit.phylum%>%
  distinct(Phylum)%>%
  tally()
# 83 families in the original data
hare_rabbit.family%>%
  distinct(Family)%>%
  tally()
# 70 genera in the original data
hare_rabbit.genus%>%
  distinct(Genus)%>%
  tally()

# 120 families in my data
ps.q.agg.family%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Family)%>%
  tally
# 241 genera in my data
ps.q.agg.genus%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Genus)%>%
  tally
# 4642 ASV in my data
ps.q.agg%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(OTU)%>%
  tally


# Check dominant taxa in my data
get_dominant_taxa_in_host(ps.q.agg.phylum,tax.rank = "Phylum",
                          host = c("hare","rabbit"),
                          nonbacterial.table = F)
# Get dominant phyla in the original dataset
hare_rabbit.phylum%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(Phylum)%>%
  summarise(mean_relab =mean(RelativeAbundance ), 
            min_relab= min(RelativeAbundance ), 
            max_relab= max(RelativeAbundance ))%>%
  arrange(desc(mean_relab)) 

# Get dominant families in the original dataset
hare_rabbit.family%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(Family)%>%
  summarise(mean_relab =mean(RelativeAbundance ), 
            min_relab= min(RelativeAbundance ), 
            max_relab= max(RelativeAbundance ))%>%
  arrange(desc(mean_relab)) %>%
  mutate(Family=sub(".*\\.\\.","",Family))%>%
  head

get_dominant_taxa_in_host(ps.q.agg.family,
                          tax.rank = "Family",
                          host = c("hare","rabbit"),
                          nonbacterial.table = F)


# Don't run the part below ############


#############
##############
# # 13. Analyse data from Debebe ####
# ### 13.1 Load data from Debebe ####
# truncationlvl<-"0" # truncation level that we chose in QIIME2
# authorname<-"biagi" # name of the folder with QIIME2 output
# read.end.type<-"single" # single reads or paired reads: decided in QIIME2
# # If you already created a workspace, just load it
# phyloseq.biagi.date_time<-"20240502_23_43_39"
# load(file.path(rdafiles.directory,paste(
#   phyloseq.biagi.date_time,
#   authorname,read.end.type,"phyloseq-summary-stats",
#   truncationlvl,
#   "workspace.RData",sep = "-")))
# 
# ### 13.1 Total number of unique ASV/phyla/families/genera per class in Debebe ####
# n.asv.per.host.biagi<-get_n_uniq_taxa_per_host(ps.q.agg.biagi,"OTU")
# n.phylum.per.host.biagi<-get_n_uniq_taxa_per_host(ps.q.agg.phylum.biagi,"Phylum")
# n.family.per.host.biagi<-get_n_uniq_taxa_per_host(ps.q.agg.family.biagi,"Family")
# n.genus.per.host.biagi<-get_n_uniq_taxa_per_host(ps.q.agg.genus.biagi,"Genus")
# 
# ### 13.2 Summary statistics table from Debebe ####
# summary.stats.table.biagi<-create_summary_stats_table(ps.q.agg.genus.biagi,
#                                                       n.asv.per.host.biagi,
#                                                       n.phylum.per.host.biagi,
#                                                       n.family.per.host.biagi,
#                                                       n.genus.per.host.biagi)
# 
# # write.table(summary.stats.table.biagi,
# #             file=file.path("./output/rtables","biagi",
# #                            paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                                  "summary-table.tsv",sep="-")),
# #             row.names = F,sep = "\t")
# 
# ### 13.3 Most abundant phyla, families, and genera in Debebe ####
# ps.q.agg.dominant.phyla.biagi<-get_dominant_taxa_in_host(ps.q.agg.phylum.biagi,
#                                                        "Phylum","NMRwt")
# ps.q.agg.dominant.families.biagi<-get_dominant_taxa_in_host(ps.q.agg.family.biagi,
#                                                           "Family","NMRwt")
# ps.q.agg.dominant.genera.nmr<-get_dominant_taxa_in_host(ps.q.agg.genus,
#                                                         "Genus","NMR")
# # save.image(paste0(rdafiles.directory,
# #                   paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                         "pooled-and-biagi-workspace.Rda")))
# 
# 
# 
# 
# ### 13.4 Bar plot of dominant families in lab vs wild NMR from our data (wild data reanalysed) ####
# ps.q.agg.dominant.families.nmr%>%
#   ungroup()%>%
#   select(Family,MeanRelativeAbundance)%>%
#   full_join(ps.q.agg.dominant.families.biagi[,c("Family","MeanRelativeAbundance")],
#             by="Family")%>%
#   rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
#   pivot_longer(cols = c("NMR","NMRwt"),
#                names_to = "class",
#                values_to = "MeanRelativeAbundance")%>%
#   arrange(-MeanRelativeAbundance)%>%
#   head(n=30)%>%
#   mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
#   ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
#              y=MeanRelativeAbundance,
#              fill=animal))+
#   geom_bar(stat = "identity",position = "dodge")+
#   theme_bw()+
#   xlab("") +
#   ylab("Average relative abundance (%)")+
#   labs(fill="")+
#   coord_cartesian(expand = FALSE)+
#   mytheme + # general ggplot theme
#   theme(axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate
#         # the x-axis labels by 45 degrees and shift to the right
#         legend.position = "inside",
#         legend.position.inside =  c(0.9,0.8)) # legend under the plot
# # ggsave(file.path("./images/barplots",
# #                  paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                        "barplot-wild-vs-lab-nmr-my-results.png")),
# #        plot=last_plot(),
# #        width = 5000,height = 4000,
# #        units = "px",dpi=300,device = "png")  
# 
# 
# 
# #### 13.4.1 Bar plot of dominant families reported in the paper ####
# wild.nmr.families<-data.frame(
#   Family=c("Lachnospiraceae",
#            "Prevotellaceae",
#            "Paraprevotellaceae",
#            "Bacteroidales Order",
#            "Clostridiales Order",
#            "Oscillospiraceae",
#            "Veillonellaceae",
#            "Clostridiaceae",
#            "Muribaculaceae",
#            "Porphyromonadaceae"),
#   MeanRelativeAbundance=c(17.6,
#                           11,
#                           8.8,
#                           6.2,
#                           6.1,
#                           5.7,
#                           4.7,
#                           4.1,
#                           4,
#                           3)
# )
# # Add a column if the family was reported as dominant or not
# ps.q.agg.dominant.families.nmr.with_domin<-ps.q.agg.dominant.families.nmr%>%
#   ungroup()%>%
#   mutate(row.index=as.numeric(rownames(ps.q.agg.dominant.families.nmr)),
#          dominant=ifelse(row.index<=10,"dominant","not_dominant"))
# 
# ps.q.agg.dominant.families.nmr.with_domin%>%
#   select(Family,MeanRelativeAbundance,dominant)%>%
#   distinct(Family,.keep_all = T)%>%
#   filter(Family%in%wild.nmr.families$Family)%>%
#   full_join(ps.q.agg.dominant.families.nmr.with_domin[1:10,c("Family","MeanRelativeAbundance","dominant")])%>%
#   full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
#             by="Family")%>%
#   rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
#   pivot_longer(cols = c("NMR","NMRwt"),
#                names_to = "class",
#                values_to = "MeanRelativeAbundance")%>%
#   arrange(-MeanRelativeAbundance)%>%
#   # head(n=30)%>%
#   mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
#   replace_na(replace =list( MeanRelativeAbundance=0))%>%
#   ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
#              y=MeanRelativeAbundance,
#              fill=animal))+
#   geom_bar(stat = "identity",position = position_dodge2())+
#   theme_bw()+
#   labs(x="",
#        y="Average relative abundance (%)")+
#   labs(fill="")+
#   coord_cartesian(expand = FALSE)+
#   mytheme + # general ggplot theme
#   theme(axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate
#         # the x-axis labels by 45 degrees and shift to the right
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.8))
# # ggsave(file.path("./images/barplots",
# #                  paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                        "barplot-wild-vs-lab-nmr-from-paper.png")),
# #        plot=last_plot(),
# #        width = 5000,height = 4000,
# #        units = "px",dpi=300,device = "png")  




