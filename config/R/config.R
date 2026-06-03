project.root <- here::here()
qiime2.results.path <-"results/1-qiime2"
# Truncation level that we chose in QIIME2:
truncationlvl<-"234"
# Name of the folder with QIIME2 output:
authorname<-"shared"
# qza_file_date_time<-"20240425_02_57_13"
qza_file_date_time<-"20260209_16_33_25"
# Single reads or paired reads (decided in QIIME2):
read.end.type<-"single"
qiime2.run_id <- paste(qza_file_date_time,read.end.type,truncationlvl,sep="-")

results<- function(...) file.path(project.root, "results","1-qiime2", 
                                  qiime2.run_id, ...)
analysis.out <- function(...) file.path(project.root, "analysis", "1-qiime2", 
                                        qiime2.run_id, ...)
maaslin2.out<- file.path(project.root, "analysis", "2-maaslin2")
manuscript.figures <- analysis.out("manuscript-figures")

analysis.params <- list(
  all_domains = list(
    label = "all_domains",
    taxa_to_keep = NULL,
    community_composition_dir = analysis.out("01-community_composition", "all_domains"),
    diversity_dir = analysis.out("02-diversity", "all_domains"),
    diffabund_dir = analysis.out("03-differential_abundance", "all_domains")
  ),
  bacteria_only = list(
    label = "bacteria",
    taxa_to_keep = "Bacteria",
    community_composition_dir = analysis.out("01-community_composition", "bacteria_only"),
    diversity_dir = analysis.out("02-diversity", "bacteria_only"),
    diffabund_dir = analysis.out("03-differential_abundance", "bacteria_only")
    
  )
)

params <- analysis.params[[active.analysis]]

community.composition.tables <-file.path(params$community_composition_dir, "tables")
community.composition.rdafiles <-file.path(params$community_composition_dir, "rdafiles")
community.composition.figures <-file.path(params$community_composition_dir, "figures")
diversity.tables <-file.path(params$diversity_dir, "tables")
diversity.rdafiles <-file.path(params$diversity_dir, "rdafiles")
diversity.figures <-file.path(params$diversity_dir, "figures")
diffabund.tables <-file.path(params$diffabund_dir, "tables")
diffabund.rdafiles <-file.path(params$diffabund_dir, "rdafiles")
diffabund.figures <-file.path(params$diffabund_dir, "figures")

# 001
qiime2.output.dir<-results("01-qiime2_output") # directory with QZA files

qiime2.metadata.dir<- ifelse(authorname == "shared", 
                             file.path("./data/metadata",
                                       authorname),
                             file.path("./data/metadata",
                                       paste(authorname,"metadata",sep = "-")))# directory with metadata

# Specify the name of your metadata file.
qiime2.metadata.filename<-file.path(qiime2.metadata.dir,
                             paste("filenames",read.end.type,
                                   authorname,"raw-supercomp.tsv", 
                                   sep = "-"))
biosample.md<-read.table(file.path(qiime2.metadata.dir,"biosample_metadata_for_ncbi.tsv"),
                         sep = "\t", header= T)

new.metadata.dir <- file.path("./data/metadata/shared")

# 002
util.functions.r <- "../utils/R"
## 2. Specifying parameters and directory/file names. #### 
# Import datasets as rds files.
ps.q.agg.asv.fname <- file.path(
  community.composition.rdafiles,
  # paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,# old
  # paste("20260211_17_01_07","phyloseq-qiime",authorname,"OTU",read.end.type, # new 
  paste("phyloseq-qiime",authorname,"OTU",read.end.type, # new 
        truncationlvl,"table.rds",sep = "-"))

ps.q.agg.genus.fname<-file.path(
  community.composition.rdafiles,
  # paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
  # paste("20260211_17_01_10","phyloseq-qiime",authorname,"Genus",read.end.type,
  paste("phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-"))

ps.q.agg.family.fname <- file.path(
  community.composition.rdafiles,
  # paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
  # paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
  paste("phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-"))

ps.q.agg.phylum.fname <- file.path(
  community.composition.rdafiles,
  # paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
  # paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
  paste("phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-"))

#' Import metadata:
custom.md.path<-file.path(new.metadata.dir,"custom.md.rds")

## Setup plots.
#' Set "pretty" labels for barplot facets that correspond to animal hosts. 
#' Here, the left side of the vector (values) is taken from the metadata, 
#' while the right side (names) are the pretty labels that will be shown on 
#' the final barplot.
pretty.level.names<-c("NMR" = "*H. glaber*", # better labels for facets
                      "DMR" = "*F. damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*N. leucodon*",
                      "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
image.formats<-c("png","tiff")



# 004
# ps.q.agg.date_time<-"20260211_17_01_10"
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts

# 005
#' Import abundance table as an rds file (NOT rarefied): 
# ps.q.agg.date_time<-"20260211_17_01_10"
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#' Import the rarefied abundance table:
# ps.q.df.pca.input.date_time<-"20260211_17_14_19"

# 006
#' Import abundance table as an rds file (NOT rarefied): 
# ps.q.agg.date_time<-"20260211_17_01_10"
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#' Import the rarefied abundance table (rds  file):
# ps.q.df.preprocessed.date_time<-"20260211_17_14_19"
# 20260211_17_14_19 rarefied table for all hosts, genus level
# 20260211_17_14_20 rarefied table file for NMR, genus level
# 20260211_17_14_21 rarefied table file for NMR, ASV level

# 008
#' Import rarefied data (rds).
#' ps.q.df.preprocessed.date_time<-"20260211_17_14_21" # ASV NMR
#' #' Import abundance table as an rds file (NOT rarefied): 
#' ps.q.agg.date_time<-"20260211_17_01_07" # ASV

# 009
#' Import datasets as rds files.
# ps.q.agg.asv<-readRDS(file=file.path(
#   rdafiles.directory,
#   # paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,# old
#   paste("20260211_17_01_07","phyloseq-qiime",authorname,"OTU",read.end.type, # new 
#         truncationlvl,"table.rds",sep = "-")))
