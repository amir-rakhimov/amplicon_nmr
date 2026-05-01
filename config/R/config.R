# 001
qiime2.results.path <-"results/1-qiime2"
#' Truncation level that we chose in QIIME2:
truncationlvl<-"234"
#' Name of the folder with QIIME2 output:
authorname<-"pooled"
# qza_file_date_time<-"20240425_02_57_13"
qza_file_date_time<-"20260209_16_33_25"
#' Single reads or paired reads (decided in QIIME2):
read.end.type<-"single"
qiime2.run_id <- paste(qza_file_date_time,read.end.type,truncationlvl,sep="-")
qiime2.output.path<-file.path(qiime2.results.path,qiime2.run_id,
                              "01-qiime2_output") # directory with QZA files

qiime2.metadata.dir<-file.path("./data/metadata",
                       paste(authorname,"metadata",sep = "-")) # directory with metadata
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)



#' Specify the name of your metadata file.
qiime2.metadata.filename<-file.path(qiime2.metadata.dir,
                             paste("filenames",read.end.type,
                                   authorname,"raw-supercomp.tsv", 
                                   sep = "-"))
biosample.md<-read.table("./data/metadata/pooled-metadata/biosample_metadata_for_ncbi.tsv",
                         sep = "\t", header= T)

# 002
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Directories with input files:

#' Import datasets as rds files.
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,# old
  paste("20260211_17_01_07","phyloseq-qiime",authorname,"OTU",read.end.type, # new 
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
  paste("20260211_17_01_10","phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
  paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
  paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

#' Import metadata:
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))

#' Specify paths and image formats:
barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")


#+ echo=FALSE
## 3. Setup plots. ####
#'
#' ## Setup plots.
#' Set "pretty" labels for barplot facets that correspond to animal hosts. 
#' Here, the left side of the vector (values) is taken from the metadata, 
#' while the right side (names) are the pretty labels that will be shown on 
#' the final barplot.
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "DMR" = "*Fukomys damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*")
#' Use only the taxa that are present in the workspace
#' (custom.md is metadata from the rdafile).
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

#' Setup general ggplot theme
mytheme<-theme(axis.text.y = element_text(size=10), # size of y axis ticks
               axis.title = element_text(size = 10), # size of axis names
               legend.text = element_text(size = 10), # size of legend text
               legend.title = element_text(size = 15), # size of legend title
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank()
)


# 003

#+ echo=FALSE
## 3. Import datasets. #### 
#'
#' ## Import datasets.
#' Import datasets as rds files.




# 004
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Path for custom metadata:
metadata.directory<-"./output/rdafiles" 
#' Import abundance table as an rds file (NOT rarefied): 
ps.q.agg.date_time<-"20260211_17_01_10"
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))

# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#' Import metadata:
custom.md<-readRDS(file.path(metadata.directory,"custom.md.rds"))
custom.levels<-c("NMR","B6mouse","MSMmouse","FVBNmouse",
                 "DMR","hare","rabbit","spalax","pvo")



# 005
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Name of the folder with QIIME2 output:
authorname<-"pooled"
#' Specify paths and image formats:
#' Import abundance table as an rds file (NOT rarefied): 
ps.q.agg.date_time<-"20260211_17_01_10"
ps.q.agg<-readRDS(file=file.path(
  "./output/rdafiles",
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#'
#' Import the rarefied abundance table:
ps.q.df.pca.input.date_time<-"20260211_17_14_19"
#' Import metadata:
custom.md<-readRDS("./output/rdafiles/custom.md.rds")
#' This is for plot filenames:
plot.types<-c("plot"="",
              "plot.with.labels"="-with-labels")
#' Set "pretty" labels
pretty.level.names<-c("NMR" = "*H. glaber*", # better labels for facets
                      "DMR" = "*F. damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*N. leucodon*",
                      "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
#' Set a general theme for ggplot2. Add more parameters depending on the plot.
mytheme <- theme(plot.title = element_text(size = 27),
                 axis.text.x = element_text(angle=0,size=20),
                 axis.text.y = element_text(size=20),
                 axis.title = element_text(size = 20),
                 legend.text = ggtext::element_markdown(size = 15),
                 legend.title = element_text(size = 25),
                 legend.position = "right",
                 plot.caption = ggtext::element_markdown(hjust = 0,
                                                         size=20),
                 plot.caption.position = "plot",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()
)


# 006
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Import abundance table as an rds file (NOT rarefied): 
ps.q.agg.date_time<-"20260211_17_01_10"
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#'
#' Import the rarefied abundance table (rds  file):
ps.q.df.preprocessed.date_time<-"20260211_17_14_19"
# 20260211_17_14_19 rarefied table for all hosts, genus level
# 20260211_17_14_20 rarefied table file for NMR, genus level
# 20260211_17_14_21 rarefied table file for NMR, ASV level

#' Import metadata:
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
#' Set "pretty" labels
pretty.level.names<-c("NMR" = "*H. glaber*", # better labels for facets
                      "DMR" = "*F. damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*N. leucodon*",
                      "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
mytheme<-theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_text(size=15),
  axis.title = element_text(size = 15),
  strip.text.x = element_text(size=15),
  plot.title = element_text(size = 15),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 18),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())

# 007
rare.status<-"rare"
filter.status<-"nonfiltered"


# 008
rare.status<-"rare"
filter.status<-"nonfiltered"
#' Import rarefied data (rds).
ps.q.df.preprocessed.date_time<-"20260211_17_14_21" # ASV NMR
#' Import abundance table as an rds file (NOT rarefied): 
ps.q.agg.date_time<-"20260211_17_01_07" # ASV
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))
#' Import metadata.
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
#' Labels for plots:
metric.labs=c('sobs'= case_when(agglom.rank == "OTU" ~ "Richness (Observed species)",
                                agglom.rank == "Genus" ~ "Richness (Observed genera)") ,
              'shannon' = "Shannon",
              'invsimpson' = "Inverse Simpson")
plot.metrics<-c("sobs","shannon",
                "invsimpson") # metrics to plot
div.indices<-c("sobs","shannon",
               "invsimpson")


# 009

#' Import datasets as rds files.
ps.q.agg.asv<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,# old
  paste("20260211_17_01_07","phyloseq-qiime",authorname,"OTU",read.end.type, # new 
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
  paste("20260211_17_01_10","phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
  paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
  paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

#' Import metadata:
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))

#' Specify paths and image formats:
barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")

#+ echo=FALSE
## 3. Setup plots. ####
#'
#' ## Setup plots.
pretty.level.names<-c("NMR" = "*H. glaber*", # better labels for facets
                      "DMR" = "*F. damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*N. leucodon*",
                      "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
