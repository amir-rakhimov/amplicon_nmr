#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 003-compare-data.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# Compare our data with results from the original papers. ####
#' # Compare our data with results from the original papers.
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' In this script, we will compare the results from original papers on DMR,
#' flying squirrels (PVO), and hares and rabbits with our 
#' re-analysis. We will also create a bar plot of dominant families in our data
#' vs the original wild NMR data.
#+ echo=FALSE
## 1. Load necessary libraries and scripts. ####
#'
#' ## Load necessary libraries and scripts.
# install.packages(c("tidyverse","ggtext"))
# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(tidyverse)
library(phyloseq)
library(ggtext)
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Name of the folder with QIIME2 output:
authorname<-"pooled" 
#' Name of the folder with QIIME2 output:
rdafiles.directory<-"./output/rdafiles"
#' Truncation level that we chose in QIIME2:
truncationlvl<-"234" 
#' Single reads or paired reads: decided in QIIME2.
read.end.type<-"single"

#+ echo=FALSE
## 3. Import datasets. #### 
#'
#' ## Import datasets.
#' Import datasets as rds files.
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(
    "20260211_17_01_07",
        "phyloseq-qiime",authorname,"OTU",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20260211_17_01_10",
        "phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
#' Import metadata:
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))

#+ echo=FALSE
## 4. Re-analysis of the DMR dataset. ####
#'
#' ## Re-analysis of the DMR dataset.
#' Data was taken from github of the original study.
dmr.asvs <- read_tsv("./data/dmr-original-data/asv_table_FreezedriedVsFrozen.tsv", 
                     col_types =cols(asv = col_character(),
                                     sample = col_character(),
                                     count = col_double(),
                                     relab = col_double() ))
dmr.metadata <- read_csv("./data/dmr-original-data/FDvsFrozenMetadata.csv", 
                         col_types = cols(.default = col_character(),
                                          SampleDate = col_date(format = ""),
                                          SampleOrder = col_number())) %>% 
  suppressWarnings()

#' Rename SampleNumber
dmr.metadata  <- dmr.metadata %>% 
  dplyr::rename(OldSampleNumber = SampleNumber, 
                SampleNumber = NewSampleNumber )
#' Make a character of the sample number
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

#' Rename unclassified taxa:
dmr.taxonomy  <- dmr.taxonomy %>%
  mutate(
    phylum = ifelse(is.na(phylum), sprintf("%s unclassified", kingdom), 
                    phylum),
    class = ifelse(is.na(class), sprintf("%s unclassified", 
                                         str_remove(phylum, "unclassified")), 
                   class),
    order = ifelse(is.na(order), sprintf("%s unclassified", 
                                         str_remove(class, "unclassified")), 
                   order),
    family = ifelse(is.na(family), sprintf("%s unclassified", 
                                           str_remove(order, "unclassified")), 
                    family),
    genus = ifelse(is.na(genus), sprintf("%s unclassified", 
                                         str_remove(family, "unclassified")), 
                   genus),
    species = ifelse(is.na(species), sprintf("%s unclassified", 
                                             str_remove(genus, "unclassified")), 
                     species))


#+ echo=FALSE
### 4.1 Filter to keep only frozen samples. ####
#'
#' ### Filter to keep only frozen samples.
#' Also, join the asv data with taxonomy table.
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
#+ echo=FALSE
### 4.2 Compare numbers of phyla, families, genera, and asv between datasets. ####
#'
#' ### Compare numbers of phyla, families, genera, and asv between datasets.
#' Phyla in original frozen (19), original total (19), and my data (21).
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Phylum")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Phylum")
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg.phylum, class=="DMR"),
                         tax.rank = "Phylum")
#' Families in original frozen (97), original total (117), and my data (156).
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Family")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Family")
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg.family, class=="DMR"),
                         tax.rank = "Family")
#' Genera in original frozen (165), original total (210), and my data (257).
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Genus")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Genus")
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg.genus, class=="DMR"),
                         tax.rank = "Genus")
#' ASVs in original frozen (1368), original total (1768), and my data (1826).
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "asv")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "asv")
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg, class=="DMR"),
                         tax.rank = "OTU")

#+ echo=FALSE
### 4.3 Classified and unclassified genera. ####
#'
#' ### Classified and unclassified genera.
#' 41 unclassified genera in frozen samples:
dmr.taxonomy.frozen[grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally
#' 51 unclassified genera in all samples:
dmr.taxonomy[grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally

#' 124 classified genera in frozen samples:
dmr.taxonomy.frozen[!grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally 
#' 159 classified genera in all samples:
dmr.taxonomy[!grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally

#+ echo=FALSE
### 4.4 How many total reads in frozen data. ####
#'
#' ### How many total reads in frozen data.
#' 1,853,295 reads in frozen data:
sum(dmr.asv.tax.df.frozen$Abundance)
#' 3,626,584 reads in total:
sum(dmr.asv.tax.df$Abundance)

#+ echo=FALSE
### 4.5 Library size in frozen data. ####
#'
#' ### Library size in frozen data.
#' Mean 92665, SD 44202
dmr.asv.tax.df.frozen%>%
  group_by(Sample)%>%
  summarise(LibrarySize=sum(Abundance))%>%
  summarise(MeanLibrarySize=mean(LibrarySize),
            SDLibrarySize=sd(LibrarySize))

#+ echo=FALSE
### 4.6 Dominant phyla in the frozen dataset vs my data (Results). ####
#'
#' ### Dominant phyla in the frozen dataset vs my data (Results).
#' Original data:
dmr.asv.tax.df.frozen %>% 
  group_by(Sample, Phylum) %>% 
  summarise(sample_relab = sum(relab)*100) %>%
  group_by(Phylum) %>% 
  summarise(MeanRelativeAbundance =mean(sample_relab), 
            MinRelativeAbundance= min(sample_relab), 
            MaxRelativeAbundance= max(sample_relab))%>%
  arrange(desc(MeanRelativeAbundance)) %>% 
  mutate(MeanRelativeAbundance = round(MeanRelativeAbundance,digits = 5),
         MinRelativeAbundance = round(MinRelativeAbundance,digits = 5),
         MaxRelativeAbundance = round(MaxRelativeAbundance,digits = 5))
#' My data:
ps.q.agg.phylum%>%
  filter(class=="DMR")%>% 
  group_by(Phylum) %>% 
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) %>% 
  mutate(MeanRelativeAbundance = round(MeanRelativeAbundance,digits = 5),
         MinRelativeAbundance = round(MinRelativeAbundance,digits = 5),
         MaxRelativeAbundance = round(MaxRelativeAbundance,digits = 5))%>%
  print(n = 21)


#+ echo=FALSE
## 5. Re-analysis of the PVO dataset. ####
#'
#' ## Re-analysis of the PVO dataset.
#' We don't have the original data, so I'm showing my results. The comparison
#' is provided in the manuscript.
#+ echo=FALSE
### 5.1 Dominant phyla (Results). ####
#'
#' ### Dominant phyla (Results).
ps.q.agg.phylum%>%
  filter(class=="pvo")%>%
  ungroup()%>%
  group_by(Phylum)%>%
  summarise(MeanRelativeAbundance = round(mean(RelativeAbundance), digits = 5), 
            MinRelativeAbundance= round(min(RelativeAbundance), digits = 5), 
            MaxRelativeAbundance= round(max(RelativeAbundance), digits = 5),
            sd_relab = round(sd(RelativeAbundance), digits =5))%>%
  arrange(desc(MeanRelativeAbundance))

#+ echo=FALSE
### 5.2 Dominant families (Results). ####
#'
#' ### Dominant families (Results).
get_dominant_taxa_in_host(ps.q.agg.family,
                          tax.rank = "Family",
                          host = "pvo",
                          nonbacterial.table = F)

#+ echo=FALSE
## 6. Re-analysis of the Hares and rabbits dataset. ####
#'
#' ## Re-analysis of the Hares and rabbits dataset.
#+ echo=FALSE
### 6.1 Remove pregnant and lactating samples. ####
#'
#' ### Remove pregnant and lactating samples.
selected.hare_rabbit.samples<-c("Hare 1","Hare 2", "Hare 3","Hare 4","Hare 5",
                                "Hare 6","Hare 7","Hare 8",
                                "Rabbit 2","Rabbit 3", "Rabbit 6","Rabbit 7",
                                "Rabbit 10", "Rabbit 11", "Rabbit 12")
#' Import abundances from the original paper (supplementary data). Only 
#' Illumina samples.
hare_rabbit.phylum<-read.table("./data/hares-rabbits-original-data/phyla.tsv",sep = "\t",
                               header = T)%>%
  as_tibble()
#' Change to long format and remove the "Unassigned" taxon.
hare_rabbit.phylum<-hare_rabbit.phylum%>%
  filter(Sample%in%selected.hare_rabbit.samples)%>%
  pivot_longer(!Sample,names_to = "Phylum",values_to = "Abundance")%>%
  filter(Phylum!="Unassigned")%>%
  filter(Abundance!=0)

#' Import family abundances.
hare_rabbit.family<-read.table("./data/hares-rabbits-original-data/families.tsv",
                               sep = "\t",
                               header = T)%>%
  as_tibble()%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Family",values_to = "Abundance")%>%
  filter(Family!="unassigned")%>%
  filter(Abundance!=0)
#' Import genus abundances. 
hare_rabbit.genus<-read.table("./data/hares-rabbits-original-data/genera.tsv",sep = "\t",
                              header = T)%>%
  as_tibble()%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Genus",values_to = "Abundance")%>%
  filter(Abundance!=0)
 
#+ echo=FALSE
### 6.2 Compare numbers of phyla, families, genera, and asv between datasets. ####
#'
#' ### Compare numbers of phyla, families, genera, and asv between datasets.
#' 14 phyla in the original data
hare_rabbit.phylum%>%
  distinct(Phylum)%>%
  tally()
#' 17 phyla in the my data
ps.q.agg.phylum%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Phylum)%>%
  tally
#' 83 families in the original data
hare_rabbit.family%>%
  distinct(Family)%>%
  tally()
#' 120 families in my data
ps.q.agg.family%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Family)%>%
  tally
#' 70 genera in the original data
hare_rabbit.genus%>%
  distinct(Genus)%>%
  tally()
#' 241 genera in my data
ps.q.agg.genus%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Genus)%>%
  tally
#' 4642 ASV in my data
ps.q.agg%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(OTU)%>%
  tally

#+ echo=FALSE
### 6.3 Compare dominant taxa between datasets. ####
#'
#' ### Compare dominant taxa between datasets.
#' Phyla in the original dataset.
hare_rabbit.phylum%>%
  mutate(class = ifelse(grepl("Hare",Sample), "hare", "rabbit"))%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(class,Phylum)%>%
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) 
#' Phyla in my dataset.
get_dominant_taxa_in_host(ps.q.agg.phylum,tax.rank = "Phylum",
                          host = c("hare","rabbit"),
                          nonbacterial.table = F)
#' Families in the original dataset.
hare_rabbit.family%>%
  mutate(class = ifelse(grepl("Hare",Sample), "hare", "rabbit"))%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(class,Family)%>%
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) %>%
  mutate(Family=sub(".*\\.\\.","",Family))%>%
  head
#' Families in my dataset.
get_dominant_taxa_in_host(ps.q.agg.family,
                          tax.rank = "Family",
                          host = c("hare","rabbit"),
                          nonbacterial.table = F)
#+ echo=FALSE
## 7. Bar plot of dominant families reported in the wild NMR paper. ####
#'
#' ## Bar plot of dominant families reported in the wild NMR paper.
wild.nmr.families<-data.frame(
  Family=c("Lachnospiraceae",
           "Prevotellaceae",
           "Paraprevotellaceae",
           "Bacteroidales Order",
           "Clostridiales Order",
           "Oscillospiraceae",
           "Veillonellaceae",
           "Clostridiaceae",
           "Muribaculaceae",
           "Porphyromonadaceae",
           "Spirochaetaceae"),
  MeanRelativeAbundance=c(17.6,
                          11,
                          8.8,
                          6.2,
                          6.1,
                          5.7,
                          4.7,
                          4.1,
                          4,
                          3,
                          10.9)
)

ps.q.agg.dominant.families.nmr<-ps.q.agg.family%>%
  filter(class=="NMR")%>%
  distinct(class,Phylum,Family, MeanRelativeAbundance)%>%
  arrange(class,desc(MeanRelativeAbundance))
  
# Add a column if the family was reported as dominant or not
ps.q.agg.dominant.families.nmr.with_domin<-ps.q.agg.dominant.families.nmr%>%
  mutate(row.index=as.numeric(rownames(ps.q.agg.dominant.families.nmr)),
         dominant=ifelse(row.index<=10,"dominant","not_dominant"))

debebe.comparison.plot<-ps.q.agg.dominant.families.nmr.with_domin%>%
  select(Family,MeanRelativeAbundance,dominant)%>%
  filter(Family%in%wild.nmr.families$Family)%>%
  full_join(ps.q.agg.dominant.families.nmr.with_domin[1:10,
                                                      c("Family",
                                                        "MeanRelativeAbundance",
                                                        "dominant")])%>%
  full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
            by="Family")%>%
  rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
  pivot_longer(cols = c("NMR","NMRwt"),
               names_to = "class",
               values_to = "MeanRelativeAbundance")%>%
  arrange(-MeanRelativeAbundance)%>%
  mutate(Family = ifelse(!grepl(" ",Family),
                         paste0("<i>",Family,"</i>"),
                         Family),
         Family = gsub ("Bacteroidales Order", "<i>Bacteroidales Order</i>", Family),
         Family = gsub ("Clostridiales Order", "<i>Clostridiales Order</i>", Family)
         )%>%
  mutate(Family=factor(Family,levels=unique(Family)))%>%
  group_by(Family)%>%
  mutate(Family_id=cur_group_id(),
         Family_id=factor(Family_id,levels=unique(sort(Family_id))))%>%
  ungroup%>%
  mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
  replace_na(replace =list( MeanRelativeAbundance=0))%>%
  mutate(Family=gsub("_", " ", Family),
         Family = stringr::str_wrap(Family,width=30),
         Family = stringr::str_replace_all(Family, "\n", "<br>"),
         Family=paste(Family_id,Family,sep = ": "),
         Family=factor(Family, levels=unique(Family)))%>%
  arrange(animal,desc(MeanRelativeAbundance))%>%
  ggplot(aes(x=Family_id,
             y=MeanRelativeAbundance,
             fill=Family))+
  geom_bar(stat = "identity",
           position = position_dodge2(),
           width=0.8 # distance between bars
           )+
  theme_bw()+
  facet_grid(~animal,
             scales="free_x",  # each species will have its own bars inside
             # facet (instead of all bars)
             space = "free_x")+
  labs(x="",
       y="Average relative abundance (%)")+
  coord_cartesian(expand = c("bottom"=FALSE))+
  scale_fill_viridis_d(option = "C")+
  theme(axis.text.x = element_text(angle=0,size=11,hjust=0.5,colour = "black"),# shift 
        # the x-axis labels to the right
        # legend.position = "inside",
        axis.title = element_text(size = 15), # size of axis names
        axis.text.y = element_text(size=11,color="black"),
        strip.text.x = ggtext::element_markdown(size = 14),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        legend.position.inside = c(0.9,0.8),
        legend.text = ggtext::element_markdown(size=13),
        legend.title = element_text(size = 16), # size of legend title
        plot.caption = element_text(size=15), # size of plot caption
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
        )
#+ fig.height = 6, fig.width = 11
print(debebe.comparison.plot)
# for(image.format in c("png","tiff")){
#   ggsave(file.path("./images/barplots",
#                  paste0(paste(format(Sys.time(),format="%Y%m%d"),
#                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                        "-barplot-wild-vs-lab-nmr-from-paper.",image.format)),
#        plot=debebe.comparison.plot,
#        width=11, height=6,units="in",
#        # width = 5700,height =2800, units = "px",
#        dpi=300,device = image.format)
# }
sessionInfo()
rm(list = ls(all=TRUE))
gc()