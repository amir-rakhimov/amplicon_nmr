library(phyloseq)
library(vegan)
library(tidyverse)
# library(pals)
library(patchwork)
library(Polychrome)

agglom.rank<-"Genus"
# truncationlvl<-"273-229"
truncationlvl<-"234"
pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "SPFmouse" = "SPF mouse,<br>B6",
                       "FukomysDamarensis" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pal" = "*Petaurista<br>alborufus lena*",
                       "pvo" = "*Pteromys<br>volans orii*",
                       "ppg" = "*Petaurista<br>philippensis grandis*",
                       "tx" = "*Trogopterus<br>xanthipes*"
                       # ,
                       # "rabbitcontrol"="*Oryctolagus cuniculus*",
                       # "harecontrol" = "*Lepus europaeus*",
                       # "ntccontrol" = "Non-treatment<br>control"
)
custom.levels=names(pretty.facet.labels)


load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
# Create a Taxon column where items are either Unclassified (previous rank) ####
# or agglom.rank (previous taxonomic rank)
# e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
agglom.rank.col<-which(colnames(ps.q.agg.rel.pooled) ==agglom.rank) 
ps.q.agg.rel.pooled$Taxon<-
  ifelse(grepl("Unclassified|Uncultured",ps.q.agg.rel.pooled[[agglom.rank.col]]), 
         ps.q.agg.rel.pooled[[agglom.rank.col]], 
         paste0(ps.q.agg.rel.pooled[[agglom.rank.col]]," (",
                ps.q.agg.rel.pooled[[agglom.rank.col+1]],")"))
# grepl finds which agglom.rank column items are uncultured or unclassified
# if true, keeps Unclassified (previous rank)
# if false (so, it's not Unclassified, actually has taxonomic classification), 
# converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)


# "Clean" column: Strip families from "Unclassified" ####
# Order the data frame by the higher taxonomic rank 
# (which is the "Clean" column)
ps.q.agg.rel.pooled$Clean<-
  gsub("Unclassified \\(|Uncultured \\(", "", ps.q.agg.rel.pooled[[agglom.rank.col+1]])
ps.q.agg.rel.pooled$Clean<-
  gsub("\\)", "", ps.q.agg.rel.pooled$Clean)

# Convert taxa with abundance<1% into Remainder ####
ps.q.agg.rel.pooled$Clean<-
  ifelse(ps.q.agg.rel.pooled$Abundance<1,
         "Remainder (Mean abundance < 1%)",
         ps.q.agg.rel.pooled$Clean)



# Separate Unclassified taxa from the rest and sort by agglom.rank-1 (higher rank) ####
# unclas is agg.rel.pooled dataset with Unclassified taxa
ps.q.agg.rel.pooled.unclas<-
  ps.q.agg.rel.pooled[grep("Unclassified|Uncultured", 
                          ps.q.agg.rel.pooled[[agglom.rank.col]]),]
# but it doesn't have Remainder taxa 
ps.q.agg.rel.pooled.unclas<-
  ps.q.agg.rel.pooled.unclas[!grepl("Remainder", ps.q.agg.rel.pooled.unclas$Clean),]

# clas is agg.rel.pooled dataset without unclassified taxa
ps.q.agg.rel.pooled.clas<-
  ps.q.agg.rel.pooled[!grepl("Unclassified|Uncultured", 
                            ps.q.agg.rel.pooled[[agglom.rank.col]]),]
# But no Remainders
ps.q.agg.rel.pooled.clas<-
  ps.q.agg.rel.pooled.clas[!grepl("Remainder", ps.q.agg.rel.pooled.clas$Clean),]


# Remainders ####
# only remainder taxa
ps.q.agg.rel.pooled.rem<-
  ps.q.agg.rel.pooled[grep("Remainder", ps.q.agg.rel.pooled$Clean),]
# Taxon.bp is for the barplot
ps.q.agg.rel.pooled.rem$Taxon<-ps.q.agg.rel.pooled.rem$Clean


## Order by the Clean column ####
ps.q.agg.rel.pooled.unclas<-
  ps.q.agg.rel.pooled.unclas[order(ps.q.agg.rel.pooled.unclas$Clean),]

ps.q.agg.rel.pooled.clas<-
  ps.q.agg.rel.pooled.clas[order(ps.q.agg.rel.pooled.clas$Clean),]


rem.row<-data.frame("Remainder (Mean abundance < 1%)",
                    "Remainder (Mean abundance < 1%)",
                    "Remainder (Mean abundance < 1%)",
                    0,
                    "Remainder (Mean abundance < 1%)",
                    "Remainder (Mean abundance < 1%)")
colnames(rem.row)<-colnames(ps.q.agg.rel.pooled.clas)

# Merge them back
ps.q.agg.rel.pooled<-rbind(ps.q.agg.rel.pooled.rem,ps.q.agg.rel.pooled.unclas,
                           ps.q.agg.rel.pooled.clas)




# Create a common legend ####
ps.q.legend.unclas<-ps.q.agg.rel.pooled.unclas
ps.q.legend.unclas<-ps.q.legend.unclas[!duplicated(ps.q.legend.unclas$Taxon),]
ps.q.legend.unclas<-ps.q.legend.unclas[order(ps.q.legend.unclas$Clean),]

ps.q.legend.clas<-ps.q.agg.rel.pooled.clas
ps.q.legend.clas<-ps.q.legend.clas[!duplicated(ps.q.legend.clas$Taxon),]
ps.q.legend.clas<-ps.q.legend.clas[order(ps.q.legend.clas$Clean),]


ps.q.legend<-rbind(rem.row,ps.q.legend.unclas,
                   ps.q.legend.clas)

# Plot the barplots ####
set.seed(1)
plot.cols<-createPalette(nrow(ps.q.legend),
                         seedcolors = rainbow(7))
# plot.cols<-createPalette(nrow(ps.q.legend),
#                          seedcolors = c("#B22222", "#0000FF","#006400", 
#                                         "#FF8C00", "#68228B"))
col.vec<-setNames(plot.cols,ps.q.legend$Taxon)


mainplot<-ps.q.agg.rel%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")"))%>% # add a column where sample names are together with sample sizes
  # filter(class=="NMR",Abundance!=0)%>%
  ggplot(aes(x=NewSample, y=Abundance,  
             fill=factor(Taxon.bp, levels=ps.q.legend$Taxon)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
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
  # guides(fill=guide_legend(ncol=1))+ # legend as one column
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
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different rodents"))+
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
        legend.position = "bottom")

ggsave(paste0("./images/barplots/",Sys.Date(),
              "-total-barplot-",truncationlvl,"-",
              agglom.rank,".tiff"),plot = mainplot,
       width = 11000,height = 4500,
       dpi = 300,units = "px",device = "tiff")  

ggsave(paste0("./images/barplots/",Sys.Date(),
              "-total-barplot-",truncationlvl,"-",
              agglom.rank,".png"),plot = mainplot,
       width = 11000,height = 5200,
       dpi = 300,units = "px",device = "png")  

# library(gridExtra)
# # TODO: arrange a grid
# ps.q.agg.rel<-ps.q.agg.rel%>%
#   left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
#   # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
#   mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
#   mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")")) # add a column where sample names are together with sample sizes
# 
# classes<-custom.levels
# classes<-c("NMR","SPFmouse","hare","rabbit")
# plots<-list()
# for (i in 1:length(classes)) {
#   # Subset data for each class
#   df_sub <- subset(ps.q.agg.rel, class == classes[i])
#   
#   # Create plot for each class
#   p <- ggplot(df_sub, aes(x=NewSample, y=Abundance,  
#                           fill=factor(Taxon.bp, levels=ps.q.legend$Taxon))) +
#     geom_bar(stat = "identity") +
#     ggtitle(paste0("Class ", classes[i]))+
#     scale_fill_manual(values = col.vec)+
#     theme(
#       axis.text.x = element_text(angle=45,hjust=1),
#       axis.line = element_blank()) +
#  
#     coord_cartesian(expand=FALSE) +
#     xlab("") +
#     ylab("Mean Relative Abundance (%)")+
#     
#     theme_bw()+
#     theme(
#       axis.text.x = element_text(angle=45,hjust=1),
#       axis.line = element_blank()) +
#     # ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats and SPF mice"))+
#     theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
#           axis.line = element_blank(), 
#           strip.text.x = ggtext::element_markdown(size = 20),
#           
#           axis.text.x = element_text(angle=45,size=20,hjust=1),
#           axis.text.y = element_text(size=20),
#           axis.title = element_text(size = 20),
#           plot.title = element_text(size = 25),
#           plot.caption = element_text(size=23))+
#     guides(fill=guide_legend(ncol=1))+
#     theme(legend.position = "right")
#   # Add plot to list of plots
#   plots[[i]] <- p
# }
# plots[[length(plots)]] <- plots[[length(plots)]] + 
#   guides(fill=guide_legend(ncol=1))+
#   theme(legend.position = "right",
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 25))
# 
# # Arrange plots in a grid
# grid.arrange(grobs = plots, ncol = 4)
# 
# plots[[2]]


for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.agg.rel%>%
    left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")"))%>%
    filter(class==custom.levels[i],Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[i])
  lvl.name<-gsub("<br>"," ", lvl.name)
  
  lvl.plot<-lvl.df%>%ggplot(aes(x=NewSample, y=Abundance,  
                      fill=factor(Taxon.bp, levels=ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))])))+
    # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
    geom_bar(stat = "identity")+ # barplot
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    facet_grid(~class, # separate species
               scales="free",  # each species will have its own bars inside facet (instead of all bars)
               space = "free", # bars will have same widths
               labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
    )+
    guides(fill=guide_legend(ncol=1))+ # legend as one column
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values = col.vec,breaks = names(col.vec))+
    xlab("") +
    ylab("Mean Relative Abundance (%)")+
    labs(fill="Taxon")+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", lvl.name))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          axis.line = element_blank(), 
          strip.text.x = ggtext::element_markdown(size = 20),
          
          axis.text.x = element_text(angle=45,size=20,hjust=1),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = ggtext::element_markdown(size = 25),
          plot.caption = element_text(size=23),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  ggsave(paste0("./images/barplots/",Sys.Date(),"-",
                custom.levels[i],"-barplot-",truncationlvl,"-",
                agglom.rank,".png"),plot = lvl.plot,
         width = 5000,height = 6000,
         dpi = 300,units = "px",device = "png")  
}

