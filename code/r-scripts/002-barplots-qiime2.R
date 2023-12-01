# Creating barplots ####

# After importing the QZA files into R and processing them with phyloseq,
# it's time to explore the taxonomic composition of our data.
# We will use the Polychrome package to create a custom palette for the 
# barplots.
library(phyloseq)
library(tidyverse)
library(Polychrome)
## Specifying parameters and directory/file names #### 
authorname<-"pooled" # name of the folder with QIIME2 output
agglom.rank<-"Genus" # this is the taxonomic rank that was used for agglomeration
truncationlvl<-"234" #  truncation level that we chose in QIIME2
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved

# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
load(paste0("./output/rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))

# Pretty labels for barplot facets that correspond to animal hosts. Here,
# the left side of the vector (values) is taken from the metadata, while
# the right side (names) are the pretty labels that will be shown on the final
# barplot
pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "B6mouse" = "B6 mouse",
                       "MSMmouse" = "MSM/Ms mouse",
                       "FVBNmouse" = "FVB/N <br>mouse",
                       "DMR" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus <br>cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pvo" = "*Pteromys volans orii*",
                       "NMRwt"="Wild *Heterocephalus glaber*"
)
# Set custom levels for the barplot. These are the animal hosts that will be used
# in barplots.
# Use only the taxa that are present in the workspace
# (custom.md is metadata from the rdafile)
custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
# Filter the phyloseq object to retain animal hosts from custom.levels
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)

## "Clean" column: Strip families from "Unclassified" ####
# Order the data frame by the higher taxonomic rank 
# (which is the "Clean" column).
# The purpose is to order our barplot legend by a higher taxonomic rank. 
# For example, if we build a barplot of genera, we may have a lot of genera
# but few families. So, for readers, it's easier to check the families first,
# and then move to genera. And when our families are ordered, it's clearly
# easier to do the checking.
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}
ps.q.agg$Clean<-
  gsub("Unclassified \\(|Uncultured \\(", "", ps.q.agg[[agglom.rank.col-1]])
ps.q.agg$Clean<-
  gsub("\\)", "", ps.q.agg$Clean)

## Convert taxa with mean relative abundance<1% into Remainder ####
# These taxa are too rare to be shown on the barplot
ps.q.agg$Clean<-
  ifelse(ps.q.agg$MeanRelativeAbundance<1,
         "Remainder (Mean abundance < 1%)",
         ps.q.agg$Clean)

## Separate Unclassified taxa from the rest and sort by agglom.rank-1 (higher rank) ####
# We want to split the dataset into three sub-datasets: unclassified taxa, 
# classified taxa, and remainders. The purpose is to order the barplots because
# they're stacked. When we split, reorder, and merge back, our barplots will have
# remainder taxa on top (for each bar), then unclassified taxa, and then finally 
# classified taxa. Moreover, we will order the data by taxonomic rank that 
# precedes the agglomerating rank (e.g. Family if we agglomerate by genera). 
# Our legend will also be ordered like the bars.

# ps.q.agg.unclas is ps.q.agg dataset with Unclassified taxa only

ps.q.agg.unclas<-
  ps.q.agg[grep("Unclassified|Uncultured", 
                          ps.q.agg[[agglom.rank.col]]),]
# But it doesn't have Remainder taxa 
ps.q.agg.unclas<-
  ps.q.agg.unclas[!grepl("Remainder", ps.q.agg.unclas$Clean),]

# clas is ps.q.agg dataset without unclassified taxa
ps.q.agg.clas<-
  ps.q.agg[!grepl("Unclassified|Uncultured", 
                            ps.q.agg[[agglom.rank.col]]),]
# But no Remainders
ps.q.agg.clas<-
  ps.q.agg.clas[!grepl("Remainder", ps.q.agg.clas$Clean),]


## Remainders ####
# Only remainder taxa
ps.q.agg.rem<-
  ps.q.agg[grep("Remainder", ps.q.agg$Clean),]
# Taxon.bp is for the barplot
ps.q.agg.rem$Taxon<-ps.q.agg.rem$Clean


### Order by the Clean column ####
ps.q.agg.unclas<-
  ps.q.agg.unclas[order(ps.q.agg.unclas$Clean),]

ps.q.agg.clas<-
  ps.q.agg.clas[order(ps.q.agg.clas$Clean),]


# Merge them back into a new ps.q.agg
ps.q.agg<-rbind(ps.q.agg.rem,ps.q.agg.unclas,
                           ps.q.agg.clas)

## Create a common legend ####
# We just copy the ps.q.agg.unclas, ps.q.agg.clas, and ps.q.agg.rem, but we remove
# duplicate taxa and order by the Clean column
ps.q.legend.unclas<-ps.q.agg.unclas
ps.q.legend.unclas<-ps.q.legend.unclas[!duplicated(ps.q.legend.unclas$Taxon),]
ps.q.legend.unclas<-ps.q.legend.unclas[order(ps.q.legend.unclas$Clean),]

ps.q.legend.clas<-ps.q.agg.clas
ps.q.legend.clas<-ps.q.legend.clas[!duplicated(ps.q.legend.clas$Taxon),]
ps.q.legend.clas<-ps.q.legend.clas[order(ps.q.legend.clas$Clean),]

# All remainder taxa will be mapped to a single row on the legend
# Here, we use `sapply` function to find remainder taxa. All text entries
# in the row such as class, Sample, or taxonomic ranks will be 
# substituted by "Remainder (Mean abundance < 1%)". Set numeric entries as zero.
ps.q.legend.rem<-ps.q.agg.clas[1,] # take the first row of ps.q.agg.clas
ps.q.legend.rem[which(sapply(ps.q.legend.rem,is.character))]<-"Remainder (Mean abundance < 1%)" # substitute all text entries with "remainder" string
ps.q.legend.rem[which(sapply(ps.q.legend.rem,is.numeric))]<-0 # set numeric
# entries as zero

# Bind sub-datasets into a new legend. For legend, select only three columns:
#   `Taxon`, `Taxon.bp`, and `Clean`
ps.q.legend<-rbind(ps.q.legend.rem,ps.q.legend.unclas,
                   ps.q.legend.clas)
ps.q.legend<-ps.q.legend%>%
  select(Taxon,Taxon.bp,Clean) # for legend

## Plot the barplots ####
# We need to choose colors for the taxa in our barplot. They should be 
# distinguishable, so we can't choose similar colors. Or at least we shouldn't 
# put them next to each other.
# 
# We will use `createPalette` function from the `Polychrome` package.
# The function takes the number of colors for the palette and seed colors that 
# will be used for palette generation. In this case, we use the number of
# rows in our legend (taxa) and rainbow colors (to avoid having similar colors 
#                                               next to each other). 
# We also need to set the random seed because the output
# is a bit random. The output is a vector of colors.

set.seed(1)
plot.cols<-createPalette(nrow(ps.q.legend),
                         seedcolors = rainbow(7))# input: number of rows
# in our legend and the seed colors that we decide to be rainbow

# The vector of colors should be named according to our legend
# because we will use this mapping in the barplot. So, each color in the vector
# correspond to each color in the legend. Remember, remainder taxa are 
# all merged into a single entry ("Remainder"), so there's just one color for 
# remainder portion.
col.vec<-setNames(plot.cols,ps.q.legend$Taxon)

# Create the barplot with ggplot2. First, we take the agglomerated
# dataset that we obtained in the `001-phyloseq-qiime2.R` and merge it with
# the dataset of total abundances, so we can know how many reads were in
# each sample. We will concatenate the sample name on the x-axis with the 
# number of reads in that sample. It will look like "`Sample 1 (n = 25000)`".
# Actually, this kind of string will be stored in a new column 
# called `NewSample`.
# 
# Then, we will convert the `class` column (host names) into factors that we 
# defined in `custom.levels`. This will order our bars according to the vector of
# levels. It must be factor because this allows us ordering `custom.levels` as
# we want, otherwise it would be alphabetic. And in our vector, the first level 
# is naked mole-rats. So, the first facet will also be naked mole-rats.
# Facets are panels that correspond to animal host.
# 
# The `ggplot` command needs an aesthetic: the x axis will correspond to
# the `NewSample` column (sample names with sample size), while the y-axis
# will be the column of relative abundances. We also need a `fill` value
# which is basically the vector that will be used for coloring the barplot.
# We use taxa from the `Taxon.bp` column because each section of each bar
# is a taxon that must be colored. But we must convert the `Taxon.bp` into
# factor, so it can map to the vector of color. **The order of factorised
# `Taxon.bp` is based on the `Taxon` column from the legend**.
mainplot<-ps.q.agg%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data 
  # with the dataset of total abundances, so we can have info about sample size
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of
  # our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ", TotalAbundance, ")"))%>% # add a 
  # column where sample names are together with sample sizes
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=ps.q.legend$Taxon)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  facet_grid(~class, # separate animal hosts
             scales="free",  # each species will have its own bars inside
             # facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels))+ # labeller will
  # change facet labels to custom
  # guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  # custom palette
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon",
       caption="Mean Relative Abundance was calculated for each host separately")+
  theme_bw()+
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different rodents"))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), #TODO: what does it do?
        strip.text.x = ggtext::element_markdown(size = 20),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "bottom") # legend under the plot
ggsave(paste0(barplot.directory,
              paste(Sys.Date(),"barplot",paste(custom.levels,collapse = '-'),
                    truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=mainplot,
       width = 13500,height = 5200,
       units = "px",dpi=300,device = "png")
ggsave(paste0(barplot.directory,
              paste(Sys.Date(),"barplot",paste(custom.levels,collapse = '-'),
                    truncationlvl,
                    agglom.rank,sep = "-"),".tiff"),
       plot=mainplot,
       width = 13500,height = 5200,
       units = "px",dpi=300,device = "tiff")

# Plot separate barplots for each host
for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.agg%>% #lvl.df is ps.q.agg. that was narrowed down
    # to the specific animal host
    left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data 
    # with the dataset of total abundances, so we can have info about sample size
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of 
    # our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%# add a 
    # column where sample names are together with sample sizes
    filter(class==custom.levels[i],Abundance!=0) # keep only rows that
  # correspond to a chosen host
  lvl.name<-unname(pretty.facet.labels[i]) # We find the pretty name for the 
  # facet using the `pretty.facet.labels` vector. `unname` will remove
  # the name from the vector element (name was taken from custom.levels, 
  # not pretty)
  lvl.name<-gsub("<br>"," ", lvl.name) # also remove all line breaks
  # the total legend is big, we need to narrow down to our host. 
  # Take the legend and extract taxa that are present in the lvl.df
  host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))]
  
  
  lvl.plot<-lvl.df%>%
    ggplot(aes(x=NewSample, y=RelativeAbundance,  
               fill=factor(Taxon.bp, levels=host.legend)))+
    # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
    geom_bar(stat = "identity")+ # barplot
    guides(fill=guide_legend(ncol=1))+ # legend as one column
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values = col.vec,
                      breaks = names(col.vec))+# custom fill that is 
    # based on our custom palette
    xlab("") +
    ylab("Relative Abundance (%)")+
    labs(fill="Taxon")+
    theme_bw()+
    ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", lvl.name))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          axis.line = element_blank(), 
          strip.text.x = ggtext::element_markdown(size = 20),# the name of 
          # each facet will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
          # the x-axis labels by 45 degrees and shift to the right
          axis.text.y = element_text(size=20), # size of y axis ticks
          axis.title = element_text(size = 20), # size of axis names
          plot.title = ggtext::element_markdown(size = 25), # the plot 
          # title will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          plot.caption = element_text(size=23),# size of plot caption
          legend.text = element_text(size = 20),# size of legend text
          legend.title = element_text(size = 25), # size of legend title
          legend.position = "right") # legend on the right
  ggsave(paste0(barplot.directory,
                paste(Sys.Date(),custom.levels[i],"barplot",truncationlvl,
                      agglom.rank,sep = "-"),".png"),
         plot=lvl.plot,
         width = 8000,height = 6000,
         units = "px",dpi=300,device = "png")
}

## Barplot for NMR and mice ####
lvl.df<-ps.q.agg%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data 
  # with the dataset of total abundances, so we can have info about sample size
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of 
  # our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
  filter(class%in%c("NMR","B6mouse"),Abundance!=0)
lvl.name<-unname(pretty.facet.labels[names(pretty.facet.labels)%in%c("NMR","B6mouse")])
lvl.name<-gsub("<br>"," ", lvl.name)
# the total legend is big, we need to narrow down to our host
host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))]

lvl.plot<-lvl.df%>%
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=host.legend)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
  )+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,breaks = names(col.vec))+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", paste(lvl.name,collapse = " and ")))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        plot.title = ggtext::element_markdown(size = 25),
        plot.caption = element_text(size=23),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")

ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-B6mouse",truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=lvl.plot,
       width = 8000,height = 6000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-B6mouse",truncationlvl,
                    agglom.rank,sep = "-"),".tiff"),
       plot=lvl.plot,
       width = 8000,height = 6000,
       units = "px",dpi=300,device = "tiff")



## Barplot for NMR and NMR wt ####
lvl.df<-ps.q.agg%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
  filter(class%in%c("NMR","NMRwt"),Abundance!=0)
lvl.name<-unname(pretty.facet.labels[names(pretty.facet.labels)%in%c("NMR","NMRwt")])
lvl.name<-gsub("<br>"," ", lvl.name)
# the total legend is big, we need to narrow down to our host
host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))]

lvl.plot<-lvl.df%>%ggplot(aes(x=NewSample, y=RelativeAbundance,  
                              fill=factor(Taxon.bp, levels=host.legend)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
  )+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,breaks = names(col.vec))+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", paste(lvl.name,collapse = " and ")))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        plot.title = ggtext::element_markdown(size = 25),
        plot.caption = element_text(size=23),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")

ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-NMRwt",truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=lvl.plot,
       width = 9000,height = 6000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-NMRwt",truncationlvl,
                    agglom.rank,sep = "-"),".tiff"),
       plot=lvl.plot,
       width = 9000,height = 6000,
       units = "px",dpi=300,device = "tiff")

# Session Info
sessionInfo()