#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
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
#' hares, and rabbits with our 
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
library(microViz)
#+ echo=FALSE
## 2. Import the config file. ####
#' 
#' ## Import the config file.
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

source(file.path(util.functions.r,"get_mean_and_sd_relab_from_phyloseq.R"))
# source(file.path(util.functions,"R/get_dominant_taxa_in_host.R"))
# source(file.path(util.functions,"R/get_n_uniq_taxa_per_host.R"))

## 3. Import datasets. #### 
#'
#' ## Import datasets.
#' Import datasets as rds files.
# ps.q.agg.asv<-readRDS(file = ps.q.agg.asv.fname)
# ps.q.agg.genus <- readRDS(file = ps.q.agg.genus.fname)
# ps.q.agg.family <- readRDS(file = ps.q.agg.family.fname)
# ps.q.agg.phylum<-readRDS(file = ps.q.agg.phylum.fname)
ps.q<-readRDS(file= ps.q.filtered.fname)
ps.q.rel <-readRDS(file= ps.q.rel.filtered.fname)
custom.md <- readRDS(custom.md.path)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

#+ echo=FALSE
## 4. Re-analysis of the DMR dataset. ####
#'
#' ## Re-analysis of the DMR dataset.
#' Data was taken from github of the original study.
dmr_original.asvs <- read_tsv("./data/dmr-original-data/asv_table_FreezedriedVsFrozen.tsv", 
                     col_types =cols(asv = col_character(),
                                     sample = col_character(),
                                     count = col_double(),
                                     relab = col_double() ))%>%
  select(-relab)%>%
  pivot_wider(names_from = "sample", 
              values_from = "count",
              values_fill = 0)%>%
  as.data.frame()%>%
  column_to_rownames("asv")
#' Convert the ASV table into phyloseq otu_table() object
dmr_original.asvs <- otu_table(dmr_original.asvs,taxa_are_rows = TRUE)

dmr_original.metadata <- read_csv("./data/dmr-original-data/FDvsFrozenMetadata.csv", 
                         col_types = cols(.default = col_character(),
                                          SampleDate = col_date(format = ""),
                                          SampleOrder = col_number())) %>% 
  suppressWarnings()%>%
  as.data.frame()%>%
  column_to_rownames("sample")
#' Convert the table into a sample_data() object
dmr_original.metadata <- sample_data(dmr_original.metadata)

#' Rename SampleNumber 
# TODO: Do we need this?
# dmr_original.metadata  <- dmr_original.metadata %>% 
#   dplyr::rename(OldSampleNumber = SampleNumber, 
#                 SampleNumber = NewSampleNumber )
#' Make a character of the sample number
# dmr_original.metadata <- dmr_original.metadata %>% 
#   mutate(Asample = paste("A", sample, sep= "_"))
dmr_original.taxonomy <- read_tsv("./data/dmr-original-data/taxonomy_FreezedriedVsFrozen.tsv", 
                         col_types =  cols(asv = col_character(),
                                           kingdom = col_character(), 
                                           phylum = col_character(),
                                           class = col_character(),
                                           order = col_character(),
                                           family = col_character(),
                                           genus = col_character(),
                                           species = col_character()))%>%
  rename("Kingdom"=kingdom)%>%
  rename("Phylum"=phylum)%>%
  rename("Class"=class)%>%
  rename("Order"=order)%>%
  rename("Family"=family)%>%
  rename("Genus"=genus)%>%
  rename("Species"=species)%>%
  as.data.frame()%>%
  column_to_rownames("asv")%>%
  as.matrix()
#' Convert the table into a tax_table() object
dmr_original.taxonomy <- tax_table(dmr_original.taxonomy)
#' Rename unclassified taxa:
dmr_original.taxonomy<-tax_fix(dmr_original.taxonomy,unknowns = c("NA","uncultured",#"Unassigned",
                                "uncultured_bacterium","uncultured_rumen",
                                "gut_metagenome","human_gut","mouse_gut",
                                "wallaby_gut","uncultured_soil", 
                                "uncultured_organism","uncultured_prokaryote",NA))

# dmr_original.taxonomy  <- dmr_original.taxonomy %>%
#   mutate(
#     phylum = ifelse(is.na(phylum), sprintf("%s unclassified", kingdom), 
#                     phylum),
#     class = ifelse(is.na(class), sprintf("%s unclassified", 
#                                          str_remove(phylum, "unclassified")), 
#                    class),
#     order = ifelse(is.na(order), sprintf("%s unclassified", 
#                                          str_remove(class, "unclassified")), 
#                    order),
#     family = ifelse(is.na(family), sprintf("%s unclassified", 
#                                            str_remove(order, "unclassified")), 
#                     family),
#     genus = ifelse(is.na(genus), sprintf("%s unclassified", 
#                                          str_remove(family, "unclassified")), 
#                    genus),
#     species = ifelse(is.na(species), sprintf("%s unclassified", 
#                                              str_remove(genus, "unclassified")), 
#                      species))

#+ echo=FALSE
### 4.1 Filter to keep only frozen samples. ####
#'
#' ### Filter to keep only frozen samples.
#' Also, join the asv data with taxonomy table.
#' Combine the data into a phyloseq object.
dmr_original.ps <- phyloseq(dmr_original.asvs,
                   dmr_original.metadata,
                   dmr_original.taxonomy)

dmr_original.ps.frozen <- subset_samples(dmr_original.ps, Treatment == "Frozen")

#' Remove taxa that have zero abundance across samples 
dmr_original.ps.frozen <- prune_taxa(taxa = taxa_names(dmr_original.ps.frozen)[(taxa_sums(dmr_original.ps.frozen) >0)],
                            x = dmr_original.ps.frozen)
dmr_original.ps <- prune_taxa(taxa = taxa_names(dmr_original.ps)[(taxa_sums(dmr_original.ps) >0)],
                            x = dmr_original.ps)
ps.q.dmr <- subset_samples(ps.q, class =="DMR")
ps.q.dmr <- prune_taxa(taxa = taxa_names(ps.q.dmr)[(taxa_sums(ps.q.dmr) >0)],
                       x = ps.q.dmr)

# dmr_original.asv.tax.df <- dmr_original.asvs%>%
#   left_join(dmr_original.taxonomy)%>%
#   rename("Sample"=sample)%>%
#   rename("Abundance"=count)%>%
#   rename("Kingdom"=kingdom)%>%
#   rename("Phylum"=phylum)%>%
#   rename("Class"=class)%>%
#   rename("Order"=order)%>%
#   rename("Family"=family)%>%
#   rename("Genus"=genus)%>%
#   rename("Species"=species)#%>%
#   # mutate(class="DMR")
# dmr_original.asv.tax.df.frozen<-dmr_original.asv.tax.df%>%
#   filter(Sample%in%subset(dmr_original.metadata,Treatment=="Frozen")$sample)
# dmr_original.taxonomy.frozen<-dmr_original.taxonomy%>%
#   filter(asv%in%dmr_original.asv.tax.df.frozen$asv)

## ====== Results after Table 1: ====== ####
#+ echo=FALSE
### 4.2 Compare numbers of phyla, families, genera, and asv between datasets. ####
#'
#' ### Compare numbers of phyla, families, genera, and asv between datasets.
#' ASVs in original frozen (1368), original total (1768), and my data (1606).
ntaxa(dmr_original.ps.frozen)
ntaxa(dmr_original.ps)
ntaxa(ps.q.dmr)
plot_bar(dmr_original.ps.frozen, fill="Phylum")
plot_bar(ps.q.dmr, fill="Phylum")

#' Phyla in original frozen (19), original total (19), and my data (20).
length(get_taxa_unique(dmr_original.ps.frozen, taxonomic.rank = "Phylum"))
length(get_taxa_unique(dmr_original.ps, taxonomic.rank = "Phylum"))
length(get_taxa_unique(ps.q.dmr, taxonomic.rank = "Phylum"))

#' Families in original frozen (97), original total (117), and my data (135).
length(get_taxa_unique(dmr_original.ps.frozen, taxonomic.rank = "Family"))
length(get_taxa_unique(dmr_original.ps, taxonomic.rank = "Family"))
length(get_taxa_unique(ps.q.dmr, taxonomic.rank = "Family"))

#' Genera in original frozen (172), original total (221), and my data (218).
length(get_taxa_unique(dmr_original.ps.frozen, taxonomic.rank = "Genus"))
length(get_taxa_unique(dmr_original.ps, taxonomic.rank = "Genus"))
length(get_taxa_unique(ps.q.dmr, taxonomic.rank = "Genus"))

#+ echo=FALSE
### 4.3 Classified and unclassified genera. ####
#'
#' ### Classified and unclassified genera.
#' 52 unclassified genera in frozen samples:
subset_taxa(dmr_original.ps.frozen, 
            grepl( "uncultured|unclassified|Kingdom|Phylum|Class|Order|Family|Ambiguous_taxa", 
                   Genus) )%>%
  get_taxa_unique(., taxonomic.rank = "Genus")%>%
  length()

#' 66 unclassified genera in all samples:
subset_taxa(dmr_original.ps, 
            grepl( "uncultured|unclassified|Kingdom|Phylum|Class|Order|Family|Ambiguous_taxa", 
                   Genus) )%>%
  get_taxa_unique(., taxonomic.rank = "Genus")%>%
  length()

#' 120 classified genera in frozen samples:
subset_taxa(dmr_original.ps.frozen, 
            !grepl( "uncultured|unclassified|Kingdom|Phylum|Class|Order|Family|Ambiguous_taxa", 
                   Genus) )%>%
  get_taxa_unique(., taxonomic.rank = "Genus")%>%
  length()
#' 155 classified genera in all samples:
subset_taxa(dmr_original.ps, 
            !grepl( "uncultured|unclassified|Kingdom|Phylum|Class|Order|Family|Ambiguous_taxa", 
                   Genus) )%>%
  get_taxa_unique(., taxonomic.rank = "Genus")%>%
  length()

#+ echo=FALSE
### 4.4 How many total reads in frozen data. ####
#'
#' ### How many total reads in frozen data.
#' 1,853,295 reads in frozen data:
sum(sample_sums(dmr_original.ps.frozen))

#' 3,626,584 reads in total:
sum(sample_sums(dmr_original.ps))

#' 2,586,498 reads in my data:
sum(sample_sums(ps.q.dmr))

#+ echo=FALSE
### 4.5 Library size in frozen data. ####
#'
#' ### Library size in frozen data.
#' Mean 92665, SD 44202
summary(sample_sums(dmr_original.ps.frozen))
sd(sample_sums(dmr_original.ps.frozen))
#' In my data, the mean is 129325, SD is 60122
summary(sample_sums(ps.q.dmr))
sd(sample_sums(ps.q.dmr))

#+ echo=FALSE
### 4.6 Dominant phyla in the frozen dataset vs my data (Results). ####
#'
#' ### Dominant phyla in the frozen dataset vs my data (Results).
#' Original data:
transform_sample_counts(dmr_original.ps.frozen, function(x) 100*x/sum(x) )%>%
  get_mean_and_sd_relab_from_phyloseq(group.name = "Treatment", 
                                      tax.rank = "Phylum")%>%
  knitr::kable(format = "simple")%>%
  print()

#' My data:
transform_sample_counts(ps.q.dmr, function(x) 100*x/sum(x) )%>%
  get_mean_and_sd_relab_from_phyloseq(group.name = "class", 
                                      tax.rank = "Phylum")%>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
## 5. Re-analysis of the Hares and rabbits dataset. ####
#'
#' ## Re-analysis of the Hares and rabbits dataset.
#+ echo=FALSE
### 5.1 Remove pregnant and lactating samples. ####
#'
#' ### Remove pregnant and lactating samples.
selected.hare_rabbit.samples<-c("Hare 1","Hare 2", "Hare 3","Hare 4","Hare 5",
                                "Hare 6","Hare 7","Hare 8",
                                "Rabbit 2","Rabbit 3", "Rabbit 6","Rabbit 7",
                                "Rabbit 10", "Rabbit 11", "Rabbit 12")

#' Filter my phyloseq object 
ps.q.hare_rabbit <- subset_samples(ps.q, class %in% c("hare","rabbit"))
ps.q.hare_rabbit <- prune_taxa(taxa = taxa_names(ps.q.hare_rabbit)[(taxa_sums(ps.q.hare_rabbit) >0)],
                       x = ps.q.hare_rabbit)

#' Import abundances from the original paper (supplementary data). Only 
#' Illumina samples.
hare_rabbit.phylum<-read.table("./data/hares-rabbits-original-data/phyla.tsv",sep = "\t",
                               header = T)%>%
  as_tibble()

#' Change to long format and remove the "Unassigned" taxon.
hare_rabbit.phylum <- hare_rabbit.phylum%>%
  filter(Sample%in%selected.hare_rabbit.samples)%>%
  pivot_longer(!Sample,names_to = "Phylum",values_to = "Abundance")%>%
  # filter(Phylum!="Unassigned")%>%
  filter(Abundance!=0)

#' Import family abundances.
hare_rabbit.family<-read.table("./data/hares-rabbits-original-data/families.tsv",
                               sep = "\t",
                               header = T)%>%
  as_tibble()%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Family",values_to = "Abundance")%>%
  # filter(Family!="unassigned")%>%
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
### 5.2 Compare numbers of phyla, families, genera, and asv between datasets. ####
#'
#' ### Compare numbers of phyla, families, genera, and asv between datasets.
#' 14 phyla in the original data
hare_rabbit.phylum%>%
  summarise(n_phylum = n_distinct(Phylum))%>%
  pull()
  
#' 14 phyla in the my data
length(get_taxa_unique(ps.q.hare_rabbit ,taxonomic.rank = "Phylum"))

#' 84 families in the original data
hare_rabbit.family%>%
  summarise(n_family = n_distinct(Family))%>%
  pull() 

#' 97 families in my data
length(get_taxa_unique(ps.q.hare_rabbit ,taxonomic.rank = "Family"))

#' 70 genera in the original data
hare_rabbit.genus%>%
  summarise(n_genus = n_distinct(Genus))%>%
  pull()

#' 211 genera in my data
length(get_taxa_unique(ps.q.hare_rabbit ,taxonomic.rank = "Genus"))

#' 5035 ASV in my data
ntaxa(ps.q.hare_rabbit)

#+ echo=FALSE
### 5.3 Compare dominant taxa between datasets. ####
#'
#' ### Compare dominant taxa between datasets.
#' Phyla in the original dataset.
hare_rabbit.phylum %>%
  mutate(class = ifelse(grepl("Hare",Sample), "hare", "rabbit"))%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(class,Phylum)%>%
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) %>%
  knitr::kable(format = "simple")%>%
  print()

#' Phyla in my dataset.
transform_sample_counts(ps.q.hare_rabbit, function(x) 100*x/sum(x) )%>%
  get_mean_and_sd_relab_from_phyloseq(group.name = "class", tax.rank = "Phylum")%>%
  arrange(desc(MeanRelativeAbundance))%>%
  knitr::kable(format = "simple")%>%
  print()

#' Families in the original dataset (showing only a few).
hare_rabbit.family%>%
  mutate(class = ifelse(grepl("Hare",Sample), "hare", "rabbit"))%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(class,Family)%>%
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  ungroup() %>%
  arrange(desc(MeanRelativeAbundance)) %>%
  mutate(Family=sub(".*\\.\\.","",Family))%>%
  head()%>%
  knitr::kable(format = "simple")%>%
  print()

#' Families in my dataset (showing only a few).
transform_sample_counts(ps.q.hare_rabbit, function(x) 100*x/sum(x) )%>%
  get_mean_and_sd_relab_from_phyloseq(group.name = "class", tax.rank = "Family")%>%
  arrange(desc(MeanRelativeAbundance))%>%
  head()%>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
## 6. Bar plot of dominant families reported in the wild NMR paper. ####
#'
#' ## Bar plot of dominant families reported in the wild NMR paper.
#' Create a dataframe with family abundances from wild NMRs (based on the 
#' results and discussion).
wild.nmr.families <- data.frame(
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

#' Subset the main phyloseq object 
ps.q.rel.nmr <- subset_samples(ps.q.rel, class == "NMR")
ps.q.rel.nmr <- prune_taxa(taxa = taxa_names(ps.q.rel.nmr)[(taxa_sums(ps.q.rel.nmr) >0)],
                               x = ps.q.rel.nmr)

#' Find average relative abundances of NMR families (in my data)
mean_sd_relab.nmr_no_groups.families <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.rel.nmr, 
                                      group.name = "class",
                                      tax.rank = "Family")

#' Add a column if the family was reported as dominant or not
mean_sd_relab.nmr_no_groups.families.dominant <- mean_sd_relab.nmr_no_groups.families %>%
  mutate(row.index = as.numeric(rownames(.)),
         dominant=ifelse(row.index <= 10,"dominant","not_dominant"))

debebe.comparison.plot <- mean_sd_relab.nmr_no_groups.families.dominant%>%
  select(Family,MeanRelativeAbundance,dominant)%>%
  filter(Family%in%wild.nmr.families$Family)%>% # intersection of my and their Families
  full_join(mean_sd_relab.nmr_no_groups.families.dominant[1:10,
                                                      c("Family",
                                                        "MeanRelativeAbundance",
                                                        "dominant")])%>% # add my dominant families that are missing
  full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
            by="Family")%>% # add the rest of their data that is missing
  rename(NMRlab=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>% 
  pivot_longer(cols = c("NMRlab","NMRwt"),
               names_to = "class",
               values_to = "MeanRelativeAbundance")%>% 
  arrange(-MeanRelativeAbundance)%>%
  mutate(Family = ifelse(!grepl(" ",Family),
                         paste0("<i>",Family,"</i>"),
                         Family),
         Family = gsub ("Bacteroidales Order", "<i>Bacteroidales Order</i>", Family),
         Family = gsub ("Clostridiales Order", "<i>Clostridiales Order</i>", Family)
         )%>% # add the <i> to make names italic
  mutate(Family=factor(Family,levels=unique(Family)))%>%
  group_by(Family)%>%
  mutate(Family_id=cur_group_id(),
         Family_id=factor(Family_id,levels=unique(sort(Family_id))))%>%
  ungroup%>%
  mutate(animal=ifelse(class=="NMRlab","Lab-bred naked mole-rats",
                       "Wild naked mole-rats"))%>%
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
  facet_grid(~animal,
             scales="free_x",  # each species will have its own bars inside
             # facet (instead of all bars)
             space = "free_x")+
  labs(x="",
       y="Average relative abundance (%)")+
  scale_fill_viridis_d(option = "C")+
  coord_cartesian(expand = c("bottom"=FALSE))+
  theme_bw()+
  project_theme 

saveRDS(debebe.comparison.plot, 
        file = file.path(manuscript.figures,"debebe.comparison.plot.rds"))

debebe.comparison.plot <- debebe.comparison.plot + 
  theme(axis.text.x = element_text(angle=0,size=7,
                                   hjust=0.5,colour = "black"),# shift 
        axis.title.x=element_blank(),
        # the x-axis labels to the right
        # legend.position = "inside",
        axis.title.y = element_text(size = 10), # size of axis names
        axis.text.y = element_text(size=10,color="black"),
        strip.text.x = ggtext::element_markdown(size = 11),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        # legend.position.inside = c(0.9,0.8),
        legend.text = ggtext::element_markdown(size=9,margin = margin(l=0.02,
                                                                      unit="cm")),
        legend.title = element_text(size = 12), # size of legend title
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.spacing.y = unit(0.01, "cm"), # distant between key text
        legend.key.spacing.x = unit(0.05, "cm"), # distant between key text
        legend.box.spacing = unit(0.004,"cm"), # space between the plot and the legend box
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
  )

#+ fig.height = 6, fig.width = 11
print(debebe.comparison.plot)

for(image.format in c("png","tiff")){
  ggsave(file = paste0("barplot-wild-vs-lab-nmr-from-paper.",image.format),
         path = community.composition.figures,
       plot = debebe.comparison.plot,
       width = 170, 
       height = 90,
       units="mm",
       dpi=300,device = image.format)
}
sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()