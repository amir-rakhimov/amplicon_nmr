library(tidyverse)
library(phyloseq)
library(DT)
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
rtables.directory<-file.path("./output/rtables",authorname)

truncationlvl<-"234" # truncation level that we chose in QIIME2
# truncationlvl<-"0" # truncation level that we chose in QIIME2

authorname<-"pooled" # name of the folder with QIIME2 output
# authorname<-"biagi" # name of the folder with QIIME2 output

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


# DMR dataset ####
dmr.asvs <- read_tsv("./FreezeDriedVSFrozen-main/data/asv_table_FreezedriedVsFrozen.tsv", col_types =cols(
  asv = col_character(),
  sample = col_character(),
  count = col_double(),
  relab = col_double() ))

dmr.metadata <- read_csv("./FreezeDriedVSFrozen-main//data/FDvsFrozenMetadata.csv", col_types = cols(
  .default = col_character(),
  SampleDate = col_date(format = ""),
  SampleOrder = col_number()
)) %>% suppressWarnings()
# Rename Sample number
dmr.metadata  <- dmr.metadata %>% dplyr::rename(OldSampleNumber = SampleNumber, SampleNumber = NewSampleNumber )
# make a character of the sample number
dmr.metadata <- dmr.metadata %>% mutate(Asample = paste("A", sample, sep= "_"))

dmr.taxonomy <- read_tsv("./FreezeDriedVSFrozen-main/data/taxonomy_FreezedriedVsFrozen.tsv", col_types =  cols(
  asv = col_character(),
  kingdom = col_character(),
  phylum = col_character(),
  class = col_character(),
  order = col_character(),
  family = col_character(),
  genus = col_character(),
  species = col_character()
))


## how many unclassified ####
### phyla ####
dmr.asvs %>% select(asv) %>% unique() %>% inner_join(dmr.taxonomy %>% filter(!is.na(phylum))) %>% summarise(n()) %>% pull() / dmr.asvs %>% select(asv) %>% unique() %>%  summarise(n()) %>% pull()
### family ####
dmr.asvs %>% select(asv) %>% unique() %>% inner_join(dmr.taxonomy %>% filter(!is.na(family)))%>% summarise(n()) %>% pull() / dmr.asvs %>% select(asv) %>% unique() %>%  summarise(n()) %>% pull()
### genus ####
dmr.asvs %>% select(asv) %>% unique() %>% inner_join(dmr.taxonomy %>% filter(!is.na(genus))) %>% summarise(n()) %>% pull() / dmr.asvs %>% select(asv) %>% unique() %>%  summarise(n()) %>% pull()

# rename unclassified
dmr.taxonomy  <- dmr.taxonomy %>%
  mutate(
    phylum = ifelse(is.na(phylum), sprintf("%s unclassified", kingdom), phylum),
    class = ifelse(is.na(class), sprintf("%s unclassified", str_remove(phylum, "unclassified")), class),
    order = ifelse(is.na(order), sprintf("%s unclassified", str_remove(class, "unclassified")), order),
    family = ifelse(is.na(family), sprintf("%s unclassified", str_remove(order, "unclassified")), family),
    genus = ifelse(is.na(genus), sprintf("%s unclassified", str_remove(family, "unclassified")), genus),
    species = ifelse(is.na(species), sprintf("%s unclassified", str_remove(genus, "unclassified")), species))


## How many asvs in frozen vs freeze dried data ####
f.asv<-dmr.asvs%>%filter(sample%in%subset(dmr.metadata,Treatment=="Frozen")$sample)%>%distinct(asv)%>%pull(asv)
fd.asv<-dmr.asvs%>%filter(sample%in%subset(dmr.metadata,Treatment=="Freeze-dried")$sample)%>%distinct(asv)%>%pull(asv)

common.asv<-intersect(f.asv,fd.asv)
all.asv<-union(f.asv,fd.asv)
length(common.asv)/length(all.asv)


## How many total reads in frozen data ####
dmr.asvs%>%filter(sample%in%subset(dmr.metadata,Treatment=="Frozen")$sample)%>%summarise(sumc=sum(count))

## mean min ma lib size in frozen and freeze-dried####
dmr.asvs %>% group_by(sample) %>% summarise(count =sum(count)) %>% summarise(min =min(count), max=max(count), mean=mean(count))
# Separate by treatment
dmr.asvs %>% left_join(dmr.metadata) %>% group_by(sample, Treatment) %>%
  summarise(count =sum(count)) %>% 
  group_by(Treatment) %>% 
  summarise(min =min(count), max=max(count), mean=mean(count),sd=sd(count))

## Filter to keep only frozen samples ####
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
dmr.asv.tax.df.frozen<-dmr.asvs%>%
  left_join(dmr.taxonomy)%>%
  rename("Sample"=sample)%>%
  filter(Sample%in%subset(dmr.metadata,Treatment=="Frozen")$sample)%>%
  rename("Abundance"=count)%>%
  rename("Kingdom"=kingdom)%>%
  rename("Phylum"=phylum)%>%
  rename("Class"=class)%>%
  rename("Order"=order)%>%
  rename("Family"=family)%>%
  rename("Genus"=genus)%>%
  rename("Species"=species)%>%
  mutate(class="DMR")
dmr.taxonomy.frozen<-dmr.taxonomy%>%filter(asv%in%dmr.asv.tax.df.frozen$asv)

# Compare numbers of phyla, families, genera, and asv between datasets
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Phylum")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Phylum")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Family")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Family")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Genus")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Genus")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "asv")
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "asv")

# Classified and unclassified genera ####
dmr.taxonomy[grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 51 unclassified genus in all
dmr.taxonomy.frozen[grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 41 unclassified genus in frozen
dmr.taxonomy[!grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 159 classified genus in all
dmr.taxonomy.frozen[!grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally# 124 classified genus in frozen

## What are the most abundant phyla on average in frozen dataset####
dmr.asvs %>% left_join(dmr.taxonomy) %>%
  filter(sample%in%subset(dmr.metadata,Treatment=="Frozen")$sample)%>%
  group_by(sample, phylum) %>% 
  summarise(relab = sum(relab)) %>%
  group_by(phylum) %>% 
  summarise(mean =mean(relab), min= min(relab), max= max(relab))%>%
  arrange(desc(mean)) %>% 
  datatable() %>%
  formatRound(columns = c("mean", "min", "max"), digits = 5)
# any ASV dominating any sample?
dmr.asvs %>% left_join(dmr.taxonomy) %>% left_join(dmr.metadata) %>%  group_by(SampleNumber, Treatment, asv, phylum) %>% 
  summarise(relab = sum(relab), .groups = "drop")  %>% arrange(desc(relab)) # sample 4 high abundance one ASV both Freeze-dried and Frozen
dmr.asvs %>% filter(asv %in% (dmr.asvs %>% left_join(dmr.taxonomy) %>% left_join(dmr.metadata) %>%  group_by(SampleNumber, Treatment, asv, phylum) %>% summarise(relab = sum(relab), .groups = "drop")  %>% 
                            arrange(desc(relab)) %>% slice(1) %>% select(asv)%>% pull())) # prevalent in 5 other samples

# Dominant phyla in the original data
get_dominant_taxa_in_host(dmr.asv.tax.df.frozen,tax.rank = "Phylum",host = "DMR",nonbacterial.table = FALSE)
dmr.asv.tax.df.frozen%>%
  group_by(Sample,Phylum)%>%
  summarise(relab=sum(relab))%>%
  group_by(Phylum)%>%
  summarise(mean =mean(relab), min= min(relab), max= max(relab)) %>% 
  arrange(desc(mean))
# Dominant genera in the original data
dmr.asv.tax.df.frozen%>%
  group_by(Sample,Genus)%>%
  summarise(relab=sum(relab))%>%
  group_by(Genus)%>%
  summarise(mean =mean(relab), min= min(relab), max= max(relab)) %>% 
  arrange(desc(mean))

# Dominant phyla in my data
get_dominant_taxa_in_host(ps.q.agg.phylum,"Phylum",c("DMR"),nonbacterial.table = F)
ps.q.agg.phylum%>%
  filter(class=="DMR")%>%
  group_by(Sample,Phylum)%>%
  summarise(relab=sum(RelativeAbundance))%>%
  group_by(Phylum)%>%
  summarise(mean =mean(relab), min= min(relab), max= max(relab)) %>% 
  arrange(desc(mean))
# Dominant genera in the original data
ps.q.agg.genus%>%
  filter(class=="DMR")%>%
  group_by(Sample,Genus)%>%
  summarise(relab=sum(RelativeAbundance))%>%
  group_by(Genus)%>%
  summarise(mean =mean(relab), min= min(relab), max= max(relab)) %>% 
  arrange(desc(mean))


# PVO dataset ####
# Dominant phyla
get_dominant_taxa_in_host(ps.q.agg.phylum,tax.rank = "Phylum",host = "pvo",nonbacterial.table = F)
get_dominant_taxa_in_host(ps.q.agg.family,tax.rank = "Family",host = "pvo",nonbacterial.table = F)
ps.q.agg.phylum%>%
  filter(class=="pvo")%>%
  ungroup()%>%
  group_by(Phylum)%>%
  summarise(minab=min(RelativeAbundance),
            maxab=max(RelativeAbundance))

# Hare and rabbits ####
custom.md%>%filter(class%in%c("hare","rabbit"))%>%arrange(class,Sample)
selected.hare_rabbit.samples<-c("Hare 1","Hare 2", "Hare 3","Hare 4","Hare 5",
                                "Hare 6","Hare 7","Hare 8",
                                "Rabbit 2","Rabbit 3", "Rabbit 6","Rabbit 7",
                                "Rabbit 10", "Rabbit 11", "Rabbit 12")
hare_rabbit.phylum<-read.table("./hares-rabbits/phyla.tsv",sep = "\t",
                               header = T)
hare_rabbit.phylum<-hare_rabbit.phylum%>%filter(Sample%in%selected.hare_rabbit.samples)

hare_rabbit.phylum<-hare_rabbit.phylum[,c("Sample",names(which(colSums(subset(hare_rabbit.phylum,select=-c(Sample,Unassigned)))!=0)))]
#####
hare_rabbit.family<-read.table("./hares-rabbits/families.tsv",sep = "\t",
                               header = T)
hare_rabbit.family<-hare_rabbit.family%>%filter(Sample_ID%in%selected.hare_rabbit.samples)

hare_rabbit.family<-hare_rabbit.family[,c("Sample_ID",names(which(colSums(subset(hare_rabbit.family,select=-c(Sample_ID,unassigned)))!=0)))]

hare_rabbit.genus<-read.table("./hares-rabbits/genera.tsv",sep = "\t",
                               header = T)
hare_rabbit.genus<-hare_rabbit.genus%>%filter(Sample_ID%in%selected.hare_rabbit.samples)

hare_rabbit.genus<-hare_rabbit.genus[,c("Sample_ID",names(which(colSums(subset(hare_rabbit.genus,select=-c(Sample_ID)))!=0)))]

ps.q.agg%>%filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(OTU)%>%
  tally

ps.q.agg.genus%>%filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Genus)%>%
  tally
ps.q.agg.family%>%filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Family)%>%
  tally
get_dominant_taxa_in_host(ps.q.agg.phylum,tax.rank = "Phylum",host = c("hare","rabbit"),nonbacterial.table = F)
get_dominant_taxa_in_host(ps.q.agg.family,tax.rank = "Family",host = c("hare","rabbit"),nonbacterial.table = F)
