library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(pals)
library(patchwork)
library(Polychrome) 

# 313-229
# 284-229
# 284-203
# 273-203
# 273-198
# 265-203
# 265-198
trulvl<-"265-198"
trunclvls<-c("313-229",
             "284-229",
             "284-203",
             "273-203",
             "273-198",
             "265-203",
             "265-198")
# qzaphyloseq.script.yasuda<-"./r-scripts/001-qza-to-phyloseq-yasuda.R"
# qzaphyloseq.script.biagi<-"./r-scripts/001-qza-to-phyloseq-biagi.R"
# source(qzaphyloseq.script.yasuda)
# source(qzaphyloseq.script.biagi)
for (trulvl in trunclvls){
  load("./rdafiles/biagi-phyloseq-workspace.RData")
  truncationlvl<-trulvl
  
  load(paste0("./rdafiles/yasuda-",truncationlvl,"-phyloseq-workspace.RData"))
  
  
  # Create a Taxon column where items are either Unclassified (previous rank) ####
  # or agglom.rank (previous taxonomic rank)
  # e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
  nmr.agg.rel.pooled$Taxon<-
    ifelse(grepl("Unclassified|Uncultured",nmr.agg.rel.pooled[[1]]), 
           nmr.agg.rel.pooled[[1]], 
           paste0(nmr.agg.rel.pooled[[1]]," (",
                  nmr.agg.rel.pooled[[2]],")"))
  
  m.agg.rel.pooled$Taxon<-
    ifelse(grepl("Unclassified|Uncultured",m.agg.rel.pooled[[1]]), 
           m.agg.rel.pooled[[1]], 
           paste0(m.agg.rel.pooled[[1]]," (",
                  m.agg.rel.pooled[[2]],")"))
  
  ps.b.agg.rel.pooled$Taxon<-
    ifelse(grepl("Unclassified|Uncultured",ps.b.agg.rel.pooled[[1]]), 
           ps.b.agg.rel.pooled[[1]], 
           paste0(ps.b.agg.rel.pooled[[1]]," (",
                  ps.b.agg.rel.pooled[[2]],")"))
  # grepl finds which agglom.rank column items are uncultured or unclassified
  # if true, keeps Unclassified (previous rank)
  # if false (so, it's not Unclassified, actually has taxonomic classification), 
  # converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)
  
  
  # "Clean" column: Strip families from "Unclassified" ####
  # Order the data frame by the higher taxonomic rank 
  # (which is the "Clean" column)
  nmr.agg.rel.pooled$Clean<-
    gsub("Unclassified \\(|Uncultured \\(", "", nmr.agg.rel.pooled[[2]])
  nmr.agg.rel.pooled$Clean<-
    gsub("\\)", "", nmr.agg.rel.pooled$Clean)
  
  # Convert taxa with abundance<1% into Remainder ####
  nmr.agg.rel.pooled$Clean<-
    ifelse(nmr.agg.rel.pooled$Abundance<1,
           "Remainder (Mean abundance < 1%)",
           nmr.agg.rel.pooled$Clean)
  
  
  # Mice 
  m.agg.rel.pooled$Clean<-
    gsub("Unclassified \\(|Uncultured \\(", "", m.agg.rel.pooled[[2]])
  m.agg.rel.pooled$Clean<-
    gsub("\\)", "", m.agg.rel.pooled$Clean)
  m.agg.rel.pooled$Clean<-
    ifelse(m.agg.rel.pooled$Abundance<1,
           "Remainder (Mean abundance < 1%)",
           m.agg.rel.pooled$Clean)
  
  # NMR Biagi
  ps.b.agg.rel.pooled$Clean<-
    gsub("Unclassified \\(|Uncultured \\(", "", ps.b.agg.rel.pooled[[2]])
  ps.b.agg.rel.pooled$Clean<-
    gsub("\\)", "", ps.b.agg.rel.pooled$Clean)
  ps.b.agg.rel.pooled$Clean<-
    ifelse(ps.b.agg.rel.pooled$Abundance<1,
           "Remainder (Mean abundance < 1%)",
           ps.b.agg.rel.pooled$Clean)
  
  # Separate Unclassified taxa from the rest and sort by agglom.rank-1 (higher rank) ####
  # unclas is agg.rel.pooled dataset with Unclassified taxa
  nmr.agg.rel.pooled.unclas<-
    nmr.agg.rel.pooled[grep("Unclassified|Uncultured", 
                            nmr.agg.rel.pooled[[1]]),]
  # but it doesn't have Remainder taxa 
  nmr.agg.rel.pooled.unclas<-
    nmr.agg.rel.pooled.unclas[!grepl("Remainder", nmr.agg.rel.pooled.unclas$Clean),]
  
  # clas is agg.rel.pooled dataset without unclassified taxa
  nmr.agg.rel.pooled.clas<-
    nmr.agg.rel.pooled[!grepl("Unclassified|Uncultured", 
                              nmr.agg.rel.pooled[[1]]),]
  # But no Remainders
  nmr.agg.rel.pooled.clas<-
    nmr.agg.rel.pooled.clas[!grepl("Remainder", nmr.agg.rel.pooled.clas$Clean),]
  
  # Mice
  m.agg.rel.pooled.unclas<-
    m.agg.rel.pooled[grep("Unclassified|Uncultured",
                          m.agg.rel.pooled[[1]]),]
  m.agg.rel.pooled.unclas<-
    m.agg.rel.pooled.unclas[!grepl("Remainder", m.agg.rel.pooled.unclas$Clean),]
  
  m.agg.rel.pooled.clas<-
    m.agg.rel.pooled[!grepl("Unclassified|Uncultured", 
                            m.agg.rel.pooled[[1]]),]
  m.agg.rel.pooled.clas<-
    m.agg.rel.pooled.clas[!grepl("Remainder", m.agg.rel.pooled.clas$Clean),]
  
  
  # Biagi NMR
  ps.b.agg.rel.pooled.unclas<-
    ps.b.agg.rel.pooled[grep("Unclassified|Uncultured",
                             ps.b.agg.rel.pooled[[1]]),]
  ps.b.agg.rel.pooled.unclas<-
    ps.b.agg.rel.pooled.unclas[!grepl("Remainder", ps.b.agg.rel.pooled.unclas$Clean),]
  
  ps.b.agg.rel.pooled.clas<-
    ps.b.agg.rel.pooled[!grepl("Unclassified|Uncultured", 
                               ps.b.agg.rel.pooled[[1]]),]
  ps.b.agg.rel.pooled.clas<-
    ps.b.agg.rel.pooled.clas[!grepl("Remainder", ps.b.agg.rel.pooled.clas$Clean),]
  
  # Remainders ####
  # only remainder taxa
  nmr.agg.rel.pooled.rem<-
    nmr.agg.rel.pooled[grep("Remainder", nmr.agg.rel.pooled$Clean),]
  # TODO: why am i doing this?
  nmr.agg.rel.pooled.rem$Taxon<-nmr.agg.rel.pooled.rem$Clean
  
  m.agg.rel.pooled.rem<-
    m.agg.rel.pooled[grep("Remainder", m.agg.rel.pooled$Clean),]
  m.agg.rel.pooled.rem$Taxon<-m.agg.rel.pooled.rem$Clean
  
  ps.b.agg.rel.pooled.rem<-
    ps.b.agg.rel.pooled[grep("Remainder", ps.b.agg.rel.pooled$Clean),]
  ps.b.agg.rel.pooled.rem$Taxon<-ps.b.agg.rel.pooled.rem$Clean
  
  ## Order by the Clean column ####
  nmr.agg.rel.pooled.unclas<- 
    nmr.agg.rel.pooled.unclas[order(nmr.agg.rel.pooled.unclas$Clean),]
  nmr.agg.rel.pooled.clas<-
    nmr.agg.rel.pooled.clas[order(nmr.agg.rel.pooled.clas$Clean),]
  
  m.agg.rel.pooled.unclas<-
    m.agg.rel.pooled.unclas[order(m.agg.rel.pooled.unclas$Clean),]
  m.agg.rel.pooled.clas<-
    m.agg.rel.pooled.clas[order(m.agg.rel.pooled.clas$Clean),]
  
  ps.b.agg.rel.pooled.unclas<-
    ps.b.agg.rel.pooled.unclas[order(ps.b.agg.rel.pooled.unclas$Clean),]
  ps.b.agg.rel.pooled.clas<-
    ps.b.agg.rel.pooled.clas[order(ps.b.agg.rel.pooled.clas$Clean),]
  
  
  # Merge them back
  nmr.agg.rel.pooled<-rbind(nmr.agg.rel.pooled.rem,nmr.agg.rel.pooled.unclas,
                            nmr.agg.rel.pooled.clas)
  
  m.agg.rel.pooled<-rbind(m.agg.rel.pooled.rem,m.agg.rel.pooled.unclas,
                          m.agg.rel.pooled.clas)
  
  ps.b.agg.rel.pooled<-rbind(ps.b.agg.rel.pooled.rem,ps.b.agg.rel.pooled.unclas,
                          ps.b.agg.rel.pooled.clas)
  
  
  ## Dummy column for x axis in ggplot ####
  host<-rep("NMR",nrow(nmr.agg.rel.pooled))
  nmr.agg.rel.pooled<-cbind(nmr.agg.rel.pooled,host)
  colnames(nmr.agg.rel.pooled)[ncol(nmr.agg.rel.pooled)]<-"host"
  
  host<-rep("Mouse",nrow(m.agg.rel.pooled))
  m.agg.rel.pooled<-cbind(m.agg.rel.pooled,host)
  colnames(m.agg.rel.pooled)[ncol(m.agg.rel.pooled)]<-"host"
  
  host<-rep("NMRwt",nrow(ps.b.agg.rel.pooled))
  ps.b.agg.rel.pooled<-cbind(ps.b.agg.rel.pooled,host)
  colnames(ps.b.agg.rel.pooled)[ncol(ps.b.agg.rel.pooled)]<-"host"
  
  
  rem.row<-data.frame("Remainder (Mean abundance < 1%)",
                      "Remainder (Mean abundance < 1%)",
                      0,"Remainder (Mean abundance < 1%)",
                      "Remainder (Mean abundance < 1%)")
  colnames(rem.row)<-colnames(nmr.agg.rel.pooled.clas)
  # Create legends for each dataset ####
  nmr.legend<-rbind(rem.row,nmr.agg.rel.pooled.unclas, nmr.agg.rel.pooled.clas)
  m.legend<-rbind(rem.row,m.agg.rel.pooled.unclas,m.agg.rel.pooled.clas)
  ps.b.legend<-rbind(rem.row,ps.b.agg.rel.pooled.unclas,ps.b.agg.rel.pooled.clas)
  
  # Create a common legend ####
  all.legends<-rbind(nmr.agg.rel.pooled.clas,
                     m.agg.rel.pooled.clas#,
                     # ps.b.agg.rel.pooled.clas
                     )
  all.legends<-all.legends[!duplicated(all.legends$Taxon),]
  all.legends<-all.legends[order(all.legends$Clean),]
  
  all.legends.unc<-rbind(nmr.agg.rel.pooled.unclas,
                         m.agg.rel.pooled.unclas
                         # ,
                         # ps.b.agg.rel.pooled.unclas
                         )
  all.legends.unc<-all.legends.unc[!duplicated(all.legends.unc$Taxon),]
  all.legends.unc<-all.legends.unc[order(all.legends.unc$Clean),]
  
  all.legends<-rbind(rem.row,all.legends.unc,
                     all.legends)
  
  # Plot the barplots ####
  set.seed(1)
  plot.cols<-createPalette(nrow(all.legends),
                           seedcolors = rainbow(7))
  
  col.vec<-setNames(plot.cols,all.legends$Taxon)
  pretty.facet.labels<-c("NMR"= "*Heterocephalus glaber*", # better labels for facets
                         "control"= "SPF mouse,<br>B6"
                         # ,
                         # "NMRwt"="*Heterocephalus glaber* wild-type"
                         )
  
  ps.total<-ps.total%>%
    bind_rows(ps.b.total)
  
  mainplot<-nmr.agg.rel%>%
    bind_rows(m.agg.rel)%>% # binds two dataframes
    # bind_rows(ps.b.agg.rel)%>%
    left_join(.,ps.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
    mutate(class=factor(class,levels=c("NMR","control")))%>% # change the order of our class column, so the NMR will be first
    
    mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")"))%>% # add a column where sample names are together with sample sizes
    ggplot(aes(x=NewSample, y=Abundance,  
               fill=factor(Taxon, levels=all.legends$Taxon)))+
    geom_bar(stat = "identity")+ # barplot
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    
    # geom_text(mapping=aes(label=paste(Sample,Abundance.y),y=-1),
    #           size=3,
    #           angle=45,
    #           hjust=1)+
  
    facet_grid(~class, # separate nmr from mice
               scales="free",  # free scales to separate bars according to species
               space = "free", # bars will have same widths
               labeller = labeller(class=pretty.facet.labels))+ # labeller will change facet labels to custom
    guides(fill=guide_legend(ncol=1))+ # legend as one column
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values = col.vec)+
    xlab("") +
    ylab("Mean Relative Abundance (%)")+
    labs(fill="Taxon",
         caption="Mean Relative Abundance was calculated for each host separately")+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats and SPF mice"))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          axis.line = element_blank(), 
          strip.text.x = ggtext::element_markdown(size = 20),
          
          axis.text.x = element_text(angle=45,size=20,hjust=1),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 25),
          plot.caption = element_text(size=23),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  
  ggsave(paste0("./images/barplots/",Sys.Date(),
                "-total-barplot-",truncationlvl,"-",
                agglom.rank,".tiff"),plot = mainplot,
         width = 11000,height = 4500,
         dpi = 300,units = "px",device = "tiff")  
  
  ggsave(paste0("./images/barplots/",Sys.Date(),
                "-total-barplot-",truncationlvl,"-",
                agglom.rank,".png"),plot = mainplot,
         width = 11000,height = 4500,
         dpi = 300,units = "px",device = "png")  
  
  to.remove <- ls()
  to.remove <- c(to.remove[!grepl("^trunclvls", to.remove)], "to.remove")
  rm(list=to.remove)
}
