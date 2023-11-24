library(qiime2R)
library(tidyverse)
library(pals)
library(ggtext)

metadata<-read_q2metadata("./data/metadata.txt")
SVs<-read_qza("./data/table-dada2-no-trunc.qza")$data
taxonomy<-
  read_qza("./data/taxonomy-dada2-no-trunc.qza")$data%>%
  parse_taxonomy()

taxasums<-summarize_taxa(SVs, taxonomy)$Family
rownames(taxasums)<-gsub("\\; NA.*$","; Unclassified",
                         rownames(taxasums))
# taxasums<-taxasums[,str_sort(colnames(taxasums),numeric = TRUE)]

taxarows<-strsplit(rownames(taxasums), "; ")
taxarows.of<-sapply(taxarows,tail,2)
taxarows.of.bind<-apply(taxarows.of, 2, paste,collapse="; ")
rownames(taxasums)<-taxarows.of.bind

taxasums.perc<-make_percent(taxasums)
taxasums.perc.ab<-sort(rowMeans(taxasums.perc),decreasing = TRUE)

threshold<-0.8
taxasums.ab.num<-length(taxasums.perc.ab[taxasums.perc.ab>threshold])

# taxa_barplot(taxasums, metadata, ntoplot = taxasums.ab.num)

plotfeats<-names(taxasums.perc.ab[1:taxasums.ab.num])

# mixedrank = function(x) order(gtools::mixedorder(x)) # eric_kernfeld
# bar<-foo%>%dplyr::arrange(Taxon,mixedrank(SampleID))



fplot<-
  make_percent(taxasums) %>%
  as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  gather(-Taxon, key="SampleID", value="Abundance") %>%
  mutate(Taxon=if_else(Taxon %in% plotfeats, Taxon,
                       paste("Taxa with relative abundance <",threshold, "%"))) %>%
  group_by(Taxon, SampleID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  ungroup() %>%
  mutate(Taxon=factor(Taxon,
                      levels=rev(c(plotfeats,
                                   paste("Taxa with relative abundance <",
                                         threshold, "%"))))) %>%
  left_join(metadata) # %>%
# arrange(Taxon,mixedrank(SampleID))

samplenum<-gsub("\\_pandaseq\\.fastq\\.gz", "",fplot$SampleID)
samplenum<-as.numeric(gsub("T", "", samplenum))
fplot<-cbind(fplot,samplenum)


png(paste0("./images/barplots/",Sys.Date(),"-taxa-bar-plot-",
           "biagi-nmr.png"), height = 800,width =1500)
bplot<-
  ggplot(fplot, aes(x=samplenum, y=Abundance, fill=Taxon)) +
  geom_bar(stat="identity") +
  # theme_q2r() +
  theme_classic(base_size=8) +
  theme(panel.border = element_rect(color="black",
                                    size=1, fill=NA)) +
  theme(axis.line = element_blank(),
        strip.background = element_blank())+
  theme(
    axis.text.x = element_text(angle=0,size=13),
    axis.text.y = element_text(size=13),
    axis.title = element_text(size = 27),
    plot.title = element_text(size = 27),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = "right") +
  coord_cartesian(expand=FALSE) +
  xlab("Sample") +
  ylab("Relative Abundance (%)")+
  ggtitle(paste0("Family level relative abundance of gut microbiota samples from naked mole-rats"))+
  scale_x_continuous(labels = as.character(1:35),
                     breaks = 1:35)+
  scale_fill_manual(values=as.vector(polychrome(taxasums.ab.num+1)))
bplot
dev.off()



average.ab<-rownames_to_column(taxasums)
colnames(average.ab)[1]<-"Taxon"
average.ab<-cbind(average.ab,
                  rowMeans(make_percent(
                    average.ab[,2:ncol(average.ab)])))
colnames(average.ab)[length(average.ab)]<-"percentage.ab"
average.ab<-cbind(average.ab,"nmr")
colnames(average.ab)[length(average.ab)]<-"Species"
average.ab<-average.ab[order(average.ab$percentage.ab,
                             decreasing = TRUE),]
sliced.average.ab<-average.ab[average.ab$percentage.ab>0.8,]


plot.main<-"Family level average relative abundance of gut microbiota samples from naked mole-rats"
png(paste0("./images/barplots/",Sys.Date(),"-avg-taxa-bar-plot-",
           "biagi-nmr.png"), height = 500,width =500)
ggplot(sliced.average.ab,
       aes(x=Species,y=percentage.ab,fill=Taxon)) +
  geom_bar(stat="identity") +
  ggtitle(plot.main)+
  theme_classic(base_size=8) +
  theme(panel.border = element_rect(color="black",
                                    size=1, fill=NA)) +
  theme(axis.line = element_blank(),
        strip.background = element_blank())+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=10),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = "right",
    plot.title=element_textbox_simple(size = 25,
                                      hjust = 1),
    plot.title.position = "plot") +
  coord_cartesian(expand=FALSE) +
  xlab("") +
  ylab("Average Relative Abundance (%)")+
  
  scale_fill_manual(values=as.vector(c(polychrome(length(sliced.average.ab)-20),stepped3(20))))
dev.off()