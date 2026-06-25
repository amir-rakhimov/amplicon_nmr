project.root <- here::here()
qiime2.results.path <-"results/1-qiime2"
# Truncation level that we chose in QIIME2:
truncationlvl<-"225_225"
# Name of the folder with QIIME2 output:
authorname<-"shared"
qza_file_date_time<-"20260521_191423"
# Single reads or paired reads (decided in QIIME2):
read.end.type<-"paired_end"
qiime2.run_id <- paste(qza_file_date_time,authorname,read.end.type,sep="-")

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
qiime2.output.dir<-results("01-qiime2_output", "merged-data",truncationlvl) # directory with QZA files

qiime2.metadata.dir<- ifelse(authorname == "shared", 
                             file.path("./data/metadata",
                                       authorname),
                             file.path("./data/metadata",
                                       paste(authorname,"metadata",sep = "-")))# directory with metadata

# Specify the name of your metadata file.
qiime2.metadata.filename<-file.path(qiime2.metadata.dir,
                             paste("filenames-paired-all.tsv"))
biosample.md<-read.table(file.path(qiime2.metadata.dir,"biosample_metadata_for_ncbi.tsv"),
                         sep = "\t", header= T)

new.metadata.dir <- file.path("./data/metadata/shared")
# Specify animal hosts to be removed from the metadata
excluded.classes <- c("pvo")

# 002
# Directory with necessary functions
util.functions.r <- "../utils/R"
# Phyloseq object
ps.q.raw.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.rds",sep="-"))
ps.q.rel.raw.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.rel.raw.rds",sep="-"))
ps.q.raw.otu_table.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.otu_table.rds",sep="-"))
ps.q.raw.tax_table.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.tax_table.rds",sep="-"))
ps.q.raw.tree.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.tree.rds",sep="-"))
rare.num_samples <- 100
# For saving phyloseq objects as dataframes
ps.q.agg.asv.fname.no_ext <- paste("phyloseq-qiime",authorname,"OTU",read.end.type, # new 
                                  truncationlvl,"table",sep = "-")
ps.q.agg.genus.fname.no_ext <-paste("phyloseq-qiime",authorname,"Genus",read.end.type,
                                    truncationlvl,"table",sep = "-")
ps.q.agg.family.fname.no_ext <- paste("phyloseq-qiime",authorname,"Family",read.end.type,
                                      truncationlvl,"table",sep = "-")
ps.q.agg.phylum.fname.no_ext <- paste("phyloseq-qiime",authorname,"Phylum",read.end.type,
                                      truncationlvl,"table",sep = "-")

ps.q.filtered.fname <- gsub("raw", "filtered", ps.q.raw.fname)
ps.q.rel.filtered.fname <- gsub("raw", "filtered", ps.q.rel.raw.fname)
ps.q.filtered.otu_table.fname <- gsub("raw", "filtered", ps.q.raw.otu_table.fname)
ps.q.filtered.tax_table.fname <- gsub("raw", "filtered", ps.q.raw.tax_table.fname)
ps.q.filtered.tree.fname <- gsub("raw", "filtered", ps.q.raw.tree.fname)


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
                      # "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
image.formats<-c("png","tiff")



