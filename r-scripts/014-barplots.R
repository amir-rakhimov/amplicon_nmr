library(phyloseq)
library(vegan)
library(tidyverse)
# library(pals)
library(patchwork)
library(Polychrome)

agglom.rank<-"Genus"
truncationlvl<-"234"
pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "SPFmouse" = "SPF mouse, B6",
                       "FukomysDamarensis" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pal" = "*Petaurista alborufus lena*",
                       "pvo" = "*Pteromys volans orii*",
                       "ppg" = "*Petaurista philippensis grandis*",
                       "tx" = "*Trogopterus xanthipes*"
                       # ,
                       # "rabbitcontrol"="*Oryctolagus cuniculus*",
                       # "harecontrol" = "*Lepus europaeus*",
                       # "ntccontrol" = "Non-treatment<br>control"
                       )
custom.levels=names(pretty.facet.labels)


load("./rdafiles/yashliu-qiime2-284-203-Genus-138-1-phyloseq-workspace.RData")

# Create a Taxon column where items are either Unclassified (previous rank) ####
# or agglom.rank (previous taxonomic rank)
# e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
agglom.rank.col<-which(colnames(ps.agg.rel.pooled) ==agglom.rank) 
ps.agg.rel.pooled$Taxon<-
  ifelse(grepl("Unclassified|Uncultured",ps.agg.rel.pooled[[agglom.rank.col]]), 
         ps.agg.rel.pooled[[agglom.rank.col]], 
         paste0(ps.agg.rel.pooled[[agglom.rank.col]]," (",
                ps.agg.rel.pooled[[agglom.rank.col+1]],")"))
# grepl finds which agglom.rank column items are uncultured or unclassified
# if true, keeps Unclassified (previous rank)
# if false (so, it's not Unclassified, actually has taxonomic classification), 
# converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)


# "Clean" column: Strip families from "Unclassified" ####
# Order the data frame by the higher taxonomic rank 
# (which is the "Clean" column)
ps.agg.rel.pooled$Clean<-
  gsub("Unclassified \\(|Uncultured \\(", "", ps.agg.rel.pooled[[agglom.rank.col+1]])
ps.agg.rel.pooled$Clean<-
  gsub("\\)", "", ps.agg.rel.pooled$Clean)

# Convert taxa with abundance<1% into Remainder ####
ps.agg.rel.pooled$Clean<-
  ifelse(ps.agg.rel.pooled$Abundance<1,
         "Remainder (Mean abundance < 1%)",
         ps.agg.rel.pooled$Clean)



# Separate Unclassified taxa from the rest and sort by agglom.rank-1 (higher rank) ####
# unclas is agg.rel.pooled dataset with Unclassified taxa
ps.agg.rel.pooled.unclas<-
  ps.agg.rel.pooled[grep("Unclassified|Uncultured", 
                          ps.agg.rel.pooled[[agglom.rank.col]]),]
# but it doesn't have Remainder taxa 
ps.agg.rel.pooled.unclas<-
  ps.agg.rel.pooled.unclas[!grepl("Remainder", ps.agg.rel.pooled.unclas$Clean),]

# clas is agg.rel.pooled dataset without unclassified taxa
ps.agg.rel.pooled.clas<-
  ps.agg.rel.pooled[!grepl("Unclassified|Uncultured", 
                            ps.agg.rel.pooled[[agglom.rank.col]]),]
# But no Remainders
ps.agg.rel.pooled.clas<-
  ps.agg.rel.pooled.clas[!grepl("Remainder", ps.agg.rel.pooled.clas$Clean),]


# Remainders ####
# only remainder taxa
ps.agg.rel.pooled.rem<-
  ps.agg.rel.pooled[grep("Remainder", ps.agg.rel.pooled$Clean),]
# TODO: why am i doing this?
ps.agg.rel.pooled.rem$Taxon<-ps.agg.rel.pooled.rem$Clean


## Order by the Clean column ####
ps.agg.rel.pooled.unclas<-
  ps.agg.rel.pooled.unclas[order(ps.agg.rel.pooled.unclas$Clean),]

ps.agg.rel.pooled.clas<-
  ps.agg.rel.pooled.clas[order(ps.agg.rel.pooled.clas$Clean),]


rem.row<-data.frame("Remainder (Mean abundance < 1%)",
                    "Remainder (Mean abundance < 1%)",
                    "Remainder (Mean abundance < 1%)",
                    0,
                    "Remainder (Mean abundance < 1%)",
                    "Remainder (Mean abundance < 1%)")
colnames(rem.row)<-colnames(ps.agg.rel.pooled.clas)

# Merge them back
ps.agg.rel.pooled<-rbind(ps.agg.rel.pooled.rem,ps.agg.rel.pooled.unclas,
                           ps.agg.rel.pooled.clas)




# Create a common legend ####
ps.legend.unclas<-ps.agg.rel.pooled.unclas
ps.legend.unclas<-ps.legend.unclas[!duplicated(ps.legend.unclas$Taxon),]
ps.legend.unclas<-ps.legend.unclas[order(ps.legend.unclas$Clean),]

ps.legend.clas<-ps.agg.rel.pooled.clas
ps.legend.clas<-ps.legend.clas[!duplicated(ps.legend.clas$Taxon),]
ps.legend.clas<-ps.legend.clas[order(ps.legend.clas$Clean),]


ps.legend<-rbind(rem.row,ps.legend.unclas,
                   ps.legend.clas)

# Plot the barplots ####
set.seed(1)
plot.cols<-createPalette(nrow(ps.legend),
                         seedcolors = rainbow(7))

col.vec<-setNames(plot.cols,ps.legend$Taxon)


library(gridExtra)
mainplot<-ps.agg.rel%>%
  left_join(.,ps.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")"))%>% # add a column where sample names are together with sample sizes
  ggplot(aes(x=NewSample, y=Abundance,  
             fill=factor(Taxon, levels=ps.legend$Taxon)))+
  geom_bar(stat = "identity")+ # barplot
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  
  # geom_text(mapping=aes(label=paste(Sample,Abundance.y),y=-1),
  #           size=3,
  #           angle=45,
  #           hjust=1)+
  
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
                )+
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

# TODO: arrange a grid
ps.agg.rel<-ps.agg.rel%>%
  left_join(.,ps.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")")) # add a column where sample names are together with sample sizes

classes<-custom.levels
classes<-c("NMR","SPFmouse","hare","rabbit")
plots<-list()
for (i in 1:length(classes)) {
  # Subset data for each class
  df_sub <- subset(ps.agg.rel, class == classes[i])
  
  # Create plot for each class
  p <- ggplot(df_sub, aes(x=NewSample, y=Abundance,  
                          fill=factor(Taxon, levels=ps.legend$Taxon))) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("Class ", classes[i]))+
    scale_fill_manual(values = col.vec)+
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
 
    coord_cartesian(expand=FALSE) +
    xlab("") +
    ylab("Mean Relative Abundance (%)")+
    
    theme_bw()+
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    # ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats and SPF mice"))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          axis.line = element_blank(), 
          strip.text.x = ggtext::element_markdown(size = 20),
          
          axis.text.x = element_text(angle=45,size=20,hjust=1),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 25),
          plot.caption = element_text(size=23))+
    theme(legend.position = "none")
  # Add plot to list of plots
  plots[[i]] <- p
}
plots[[length(plots)]] <- plots[[length(plots)]] + 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

# Arrange plots in a grid
grid.arrange(grobs = plots, ncol = 3)
