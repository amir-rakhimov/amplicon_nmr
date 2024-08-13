# Creating barplots ####

# After importing the QZA files into R and processing them with phyloseq,
# it's time to explore the taxonomic composition of our data.
# We will use the Polychrome package to create a custom palette for the 
# barplots.
# install.packages(c("tidyverse","ggtext","Polychrome"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# library(phyloseq)
library(tidyverse)
library(Polychrome)
library(ggtext)
## Specifying parameters and directory/file names #### 
authorname<-"pooled" # name of the folder with QIIME2 output
agglom.rank<-"OTU" # this is the taxonomic rank that was used for agglomeration
truncationlvl<-"234" #  truncation level that we chose in QIIME2
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved
image.formats<-c("png","tiff")
# phyloseq.workspace.date_time is a workspace with ps.q.agg dataframe and metadata
# from 001-phyloseq-qiime2.R
phyloseq.workspace.date_time="20240524_13_54_21" 
# 20240426_21_44_30 for all hosts, genus level
# 20240524_13_54_21 for all hosts, OTU level
# 20240426_21_43_29 for all hosts, family level
# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
load(file.path("./output/rdafiles",paste(
  phyloseq.workspace.date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

# Pretty labels for barplot facets that correspond to animal hosts. Here,
# the left side of the vector (values) is taken from the metadata, while
# the right side (names) are the pretty labels that will be shown on the final
# barplot
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
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
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
# Filter the phyloseq object to retain animal hosts from custom.levels
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)

# Check the total number of taxa (not filtered by mean relative abundance)
ps.q.agg%>%
  ungroup()%>%
  select(matches(paste0("^",agglom.rank,"$")))%>%
  pull(.)%>%
  unique()%>%
  sort()%>%
  length()
# select(matches(paste0("^",agglom.rank,"$"))) will select variables that 
# match a pattern (regex).
# So, if we are agglomerating by Phylum, the command will select the Phylum column

# Creating a legend ####
### 1. Extract the unique taxa (corresponding to agglom.rank) from the dataset ####
# Result: taxa.list vector

taxa.list<-ps.q.agg%>%
  group_by_at(c("class",agglom.rank))%>%
  filter(MeanRelativeAbundance>=1)%>%
  ungroup()%>%
  select(matches(paste0("^",agglom.rank,"$")))%>%
  pull(.)%>%
  unique()

# select(matches(paste0("^",agglom.rank,"$"))) will select variables that 
# match a pattern (regex).
# So, if we are agglomerating by Phylum, the command will select the Phylum column

### 2. Get the taxonomic ranks for ordering ####
# The order of the legend is: "Remainder","Kingdom", "Phylum", "Class", 
# "Order", "Family". But that is when we agglomerate by Genus.
# If we agglomerate by a higher rank, we need to remove the unnecessary ranks.
# For example, if agglom.rank="Phylum", we need to match the agglom.rank
# with the custom.order vector and extract the ranks before agglom.rank
# (excluding the agglom.rank!). The resulting custom.order vector becomes shorter
taxa.list<-c(taxa.list,"Remainder")
custom_order <- c("Remainder","Kingdom", "Phylum", "Class", "Order", "Family")
# If we agglomerate by higher level (Order,Class, etc), need to adjust the rank
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}

### 3. Find classified and unclassified taxa indices. ####
# The sub command removes everything before a space (e.g. Bacteria Kingdom becomes
# Bacteria) in the taxa.list vector. The match commmand matches the taxa.list
# to custom.order (taxonomic ranks before agglom.rank). This will find 
# unclassified taxa because they are from higher taxonomic rank (e.g. if we 
# agglomerate by Genus, taxa with "Family" or "Order" in their name are 
# unclassified ones). The rank vector has length equal to length of the
# taxa.list. Each value corresponds to an element from taxa.list that matches
# to the index of custom.order vector. So, if the 2nd element of taxa.list is
# "Kingdom Bacteria", the 2nd element of the rank vector is "1" because the 
# taxonomic rank "Kingdom" matches "Kingdom Bacteria" and is 1st in the custom.order
# vector.
# However, the classified taxa won't match the custom.order 
# vector (because they are classified). Thus, the classified elements will be NA.
# We substitute NA values with the length of custom.order +1. So, if we agglomerate
# by Genus, the length of custom.order is 6. But the classified taxa (which were
# not found in custom.order and thus have NA in the rank vector) are out of the 
# ranking. Their index is 6+1=7. So, in the rank vector all NAs will become 7.
rank <- match(sub(".* ", "", taxa.list),custom_order)
rank[is.na(rank)] <- length(custom_order) + 1

### 4. Extract classified taxa: their rank is the last in the custom_order vector ####
# We match taxa.list to the rank vector and extract elements with rank values
# equal to length(custom.order)+1.
classified.taxa<-taxa.list[rank==length(custom_order) + 1]
### 4.1.Unclassified taxa are ones with indices found in the custom.order ####
unclassified.taxa<-taxa.list[rank!=length(custom_order) + 1]

### 4.3. In the ps.q.agg dataset, we find the index of the agglom.rank column and  ####
# preceding.rank column (e.g Genus and Family)
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}

### 4.4. Create a dataframe of classified taxa (agglom.rank + preceding rank) for the legend ####
# From the ps.q.agg dataset, we obtain classified taxa, select columns
# for agglom.rank and preceding.rank, then retain only unique rows
# To create the Taxon.bp column, we concatenate agglom.rank with preceding.rank
# The taxa.for_bp.df has columns for agglom.rank (e.g Genus), preceding.rank 
# (e.g. Family), and Taxon.bp (merged agglom.rank and preceding.rank)
taxa.for_bp.df<-ps.q.agg%>%
  ungroup()%>%
  filter(get(agglom.rank)%in%classified.taxa)%>%
  select(all_of(c(agglom.rank,preceding.rank)))%>%
  distinct()%>%
  unite("Taxon.bp",agglom.rank:preceding.rank,sep = " (",remove = FALSE)%>%
  select(all_of(c(agglom.rank,preceding.rank,"Taxon.bp")))%>%
  mutate("Append"=")")%>%
  unite("Taxon.bp",Taxon.bp:Append,sep = "")

# taxa.for_bp.df$Taxon[is.na(taxa.for_bp.df$Taxon)]<-taxa.for_bp.df$Genus[is.na(taxa.for_bp.df$Taxon)]
# taxa.for_bp.df$Family[is.na(taxa.for_bp.df$Family)]<-taxa.for_bp.df$Genus[is.na(taxa.for_bp.df$Family)]
### 4.5. Order classified taxa by preceding.rank then agglom.rank ####
# These are agglom.rank (preceding.rank) format strings for barplot
taxa.for_bp.df<-taxa.for_bp.df%>%
  arrange(get(preceding.rank),get(agglom.rank))
### 4.6 We select the Taxon.bp column of the taxa.for_bp.df and store it as a  ####
# separate vector taxa.for_bp.list
taxa.for_bp.list<-taxa.for_bp.df$Taxon.bp

### 5. Now we order unclassified taxa: start by matching ####
### 5.1 Match the unclassified taxa to indices of custom.order ####
newrank <- match(sub(".* ", "", unclassified.taxa),custom_order)
# newrank[is.na(newrank)] <- length(custom_order) + 1

### 5.2 Order the unclassified taxa vector ####
unclassified.taxa<-unclassified.taxa[order(newrank)]

### 5.3 Sorting inside each rank: split the vector by taxonomic rank ####
# If agglom.rank="Genus", unclassified.taxa.split will be a list of 6 lists
# (Remainder, Kingdom, Phylum, Class, "Order", "Family").
unclassified.taxa.split <- split(unclassified.taxa, newrank[order(newrank)])
# Sort inside each rank
unclassified.taxa.sorted <- unlist(lapply(unclassified.taxa.split, sort))
# The elements in unclassified.taxa.sorted are named by their rank+position inside
# rank (e.g. the 19th taxon in Family will be 619: 6 for the Family list, 19 for the 
# index). Get rid of these names!
unclassified.taxa.sorted<-unname(unclassified.taxa.sorted)

### 6. Add sorted unclassified taxa to the sorted classified taxa (taxa.for_bp.list) ####
taxa.for_bp.list<-c(unclassified.taxa.sorted,taxa.for_bp.list)
# Remainder is the first element
taxa.for_bp.list[1]<-"Remainder (Mean relative abundance < 1%)"

### 7. Add barplot taxa to the main dataframe ####
# Join the ps.q.agg dataset with classified taxa dataset (taxa.for_bp.df)
ps.q.agg.for_bp<-ps.q.agg%>%
  left_join(taxa.for_bp.df,by=agglom.rank)%>%
  ungroup()%>%
  select(-paste0(preceding.rank,".y"))%>% # remove the preceding.rank column from taxa.for_bp.df
  rename(!!preceding.rank:=paste0(preceding.rank,".x")) # rename the preceding.rank column
# !!preceding.rank:= will evaluate the variable

### 7.1 Remove NAs: unclassified taxa with MeanRelativeAbundance>=1 ####
# They are NA because they were not found in taxa.for_bp.df (which has only
# classified data)
# We remove NA by simply assigning agglom.rank to Taxon.bp. So, if Genus is 
# "Bacteria Kingdom", the Taxon.bp will also become "Bacteria Kingdom".
ps.q.agg.for_bp[which(ps.q.agg.for_bp$MeanRelativeAbundance>=1&is.na(ps.q.agg.for_bp$Taxon.bp)),"Taxon.bp"]<-
  ps.q.agg.for_bp[which(ps.q.agg.for_bp$MeanRelativeAbundance>=1&is.na(ps.q.agg.for_bp$Taxon.bp)),agglom.rank]
# Taxa with MeanRelativeAbundance<1% become "Remainder"
ps.q.agg.for_bp[which(ps.q.agg.for_bp$MeanRelativeAbundance<1),"Taxon.bp"]<-
  "Remainder (Mean relative abundance < 1%)"




### 8. We want to highlight NMR-specific taxa on the barplot #### 
### 8.1 Find all unique taxa in NMR samples ####
nmr.set<-ps.q.agg.for_bp%>%
  filter(class=="NMR")%>%
  select(all_of(agglom.rank))%>%
  unique()%>%
  pull()
### 8.2 Find all unique taxa in non-NMR samples ####
others.set<-ps.q.agg.for_bp%>%
  filter(class!="NMR")%>%
  select(all_of(agglom.rank))%>%
  unique()%>%
  pull()
### 8.3  Get NMR-specific taxa ####
nmr.uniq<-setdiff(nmr.set,others.set)
others.uniq<-setdiff(others.set,nmr.set)
### 9. Obtain taxa found in the taxa.for_bp.list vector (sorted unclassified and  ####
# classified taxa) with MeanRelativeAbundance>=1%. Keep only unique values in the
# agglom.rank.vec vector
agglom.rank.vec<-ps.q.agg.for_bp%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=1)%>%
  select(agglom.rank)%>%
  pull()%>%
  unique()

### 9.1 Find NMR-specific taxa in the agglom.rank.vec vector  ####
nmr.uniq.legend<-agglom.rank.vec[agglom.rank.vec%in%nmr.uniq]
nmr.uniq.legend<-ps.q.agg.for_bp%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=1)%>%
  distinct(get(agglom.rank),Taxon.bp)%>%
  rename(!!agglom.rank:="get(agglom.rank)")%>%
  filter(get(agglom.rank)%in%nmr.uniq)%>%
  select(Taxon.bp)%>%
  pull()
others.uniq.legend<-agglom.rank.vec[agglom.rank.vec%in%others.uniq]
others.uniq.legend<-ps.q.agg.for_bp%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=1)%>%
  distinct(get(agglom.rank),Taxon.bp)%>%
  rename(!!agglom.rank:="get(agglom.rank)")%>%
  filter(get(agglom.rank)%in%others.uniq)%>%
  select(Taxon.bp)%>%
  pull()
### 10. New font colors: NMR-specific taxa will be red in the legend
ps.q.legend<-as.data.frame(taxa.for_bp.list)%>%
  rename("Taxon.bp"="taxa.for_bp.list")%>%
  mutate(new.colors=ifelse(Taxon.bp%in%nmr.uniq.legend,
                           paste("<span style='color: red'><b>",Taxon.bp,"</b></span>"),
                           Taxon.bp),
         new.colors=ifelse(Taxon.bp%in%others.uniq.legend,
                           paste("<span style='color: blue'><b>",Taxon.bp,"</b></span>"),
                           new.colors))
## 11. Plot the barplots ####
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
plot.cols<-createPalette(length(taxa.for_bp.list),
                         seedcolors =rainbow(7))# input: number of rows

# in our legend and the seed colors that we decide to be rainbow

# The vector of colors should be named according to our legend
# because we will use this mapping in the barplot. So, each color in the vector
# correspond to each color in the legend. Remember, remainder taxa are 
# all merged into a single entry ("Remainder"), so there's just one color for 
# remainder portion.
plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,ps.q.legend$Taxon.bp)
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
mainplot<-ps.q.agg.for_bp%>%
  group_by(Sample)%>%
  mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
  # so we can have info about sample size
  ungroup()%>%
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of
  # our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ", TotalAbundance, ")"))%>% # add a 
  # column where sample names are together with sample sizes
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=ps.q.legend$Taxon.bp)))+
  geom_bar(stat = "identity")+ # barplot
  facet_grid(~class, # separate animal hosts
             scales="free",  # each species will have its own bars inside
             # facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.level.names))+ # labeller will
  # change facet labels to custom
  # guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,
                    labels=ps.q.legend$new.colors)+ # custom fill that is based on our 
  # custom palette
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon",
       caption="Mean Relative Abundance was calculated for each host separately")+
  theme_bw()+
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different rodents"))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
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
        legend.text = element_markdown(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "bottom") # legend under the plot
# ggsave(paste0(barplot.directory,
#               paste(paste(format(Sys.time(),format="%Y%m%d"),
#                           format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                     "barplot",paste(custom.levels,collapse = '-'),
#                     truncationlvl,
#                     agglom.rank,sep = "-"),".png"),
#        plot=mainplot,
#        width = 13500,height = 5200,
#        units = "px",dpi=300,device = "png")
# ggsave(paste0(barplot.directory,
#               paste(paste(format(Sys.time(),format="%Y%m%d"),
#                           format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                     "barplot",paste(custom.levels,collapse = '-'),
#                     truncationlvl,
#                     agglom.rank,sep = "-"),".tiff"),
#        plot=mainplot,
#        width = 13500,height = 5200,
#        units = "px",dpi=300,device = "tiff")
for(image.format in image.formats){
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "barplot",paste(custom.levels,collapse = '-'),
                      truncationlvl,agglom.rank,
                      sep = "-"),".",image.format),
         plot=mainplot,
         width = 13500,height = 5200,
         units = "px",dpi=300,device = image.format)
}
# Plot separate barplots for each host ####
for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.agg.for_bp%>% #lvl.df is ps.q.agg.for_bp. that was narrowed down
    # to the specific animal host
    group_by(Sample)%>%
    mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
    # so we can have info about sample size
    ungroup()%>%
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of 
    # our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%# add a 
    # column where sample names are together with sample sizes
    filter(class==custom.levels[i],Abundance!=0) # keep only rows that
  # correspond to a chosen host
  lvl.name<-unname(pretty.level.names[i]) # We find the pretty name for the 
  # facet using the `pretty.level.names` vector. `unname` will remove
  # the name from the vector element (name was taken from custom.levels, 
  # not pretty)
  lvl.name<-gsub("<br>"," ", lvl.name) # also remove all line breaks
  # the total legend is big, we need to narrow down to our host. 
  # Take the legend and extract taxa that are present in the lvl.df
  host.legend<-ps.q.legend$Taxon.bp[ps.q.legend$Taxon.bp%in%names(table(lvl.df$Taxon.bp))]
  
  
  lvl.plot<-lvl.df%>%
    ggplot(aes(x=NewSample, y=RelativeAbundance,  
               fill=factor(Taxon.bp, levels=host.legend)))+
    geom_bar(stat = "identity")+ # barplot
    guides(fill=guide_legend(ncol=1))+ # legend as one column
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values = col.vec[names(col.vec)%in%host.legend],
                      breaks = names(col.vec)[names(col.vec)%in%host.legend],
                      labels=host.legend)+# custom fill that is 
    # based on our custom palette
    xlab("") +
    ylab("Relative Abundance (%)")+
    labs(fill="Taxon")+
    theme_bw()+
    ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", lvl.name))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
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
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      custom.levels[i],"barplot",truncationlvl,
                      agglom.rank,sep = "-"),".png"),
         plot=lvl.plot,
         width = 8000,height = 6000,
         units = "px",dpi=300,device = "png")
}

## Barplot for NMR and mice ####
lvl.df<-ps.q.agg.for_bp%>%
  group_by(Sample)%>%
  mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
  # so we can have info about sample size
  ungroup()%>%
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of 
  # our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
  filter(class%in%c("NMR","B6mouse"),Abundance!=0)
lvl.name<-unname(pretty.level.names[names(pretty.level.names)%in%c("NMR","B6mouse")])
lvl.name<-gsub("<br>"," ", lvl.name)
# the total legend is big, we need to narrow down to our host
host.legend<-ps.q.legend$Taxon.bp[ps.q.legend$Taxon.bp%in%names(table(lvl.df$Taxon.bp))]

nmr.set<-lvl.df%>%
  filter(class=="NMR")%>%
  select(agglom.rank)%>%
  unique()%>%
  pull()

b6.set<-lvl.df%>%
  filter(class=="B6mouse")%>%
  select(agglom.rank)%>%
  unique()%>%
  pull()

nmr.uniq<-setdiff(nmr.set,b6.set)
nmr.uniq.legend<-ps.q.agg.for_bp%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=1)%>%
  distinct(get(agglom.rank),Taxon.bp)%>%
  rename(!!agglom.rank:="get(agglom.rank)")%>%
  filter(get(agglom.rank)%in%nmr.uniq)%>%
  select(Taxon.bp)%>%
  pull()

b6.uniq<-setdiff(b6.set,nmr.set)
b6.uniq.legend<-ps.q.agg.for_bp%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=1)%>%
  distinct(get(agglom.rank),Taxon.bp)%>%
  rename(!!agglom.rank:="get(agglom.rank)")%>%
  filter(get(agglom.rank)%in%b6.uniq)%>%
  select(Taxon.bp)%>%
  pull()


# New font colors
host.legend<-data.frame(host.legend,host.legend)%>%
  rename(old.colors="host.legend",
         new.colors="host.legend.1")%>%
  mutate(new.colors=ifelse(host.legend%in%nmr.uniq.legend,
                           paste("<span style='color: red'><b>",new.colors,"</b></span>"),
                           old.colors),
         new.colors=ifelse(host.legend%in%b6.uniq.legend,
                           paste("<span style='color: blue'><b>",new.colors,"</b></span>"),
                           new.colors))

lvl.plot<-lvl.df%>%
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=host.legend$old.colors)))+
  geom_bar(stat = "identity")+ # barplot
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.level.names) # labeller will change facet labels to custom
  )+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec[names(col.vec)%in%host.legend$old.colors],
                    breaks = names(col.vec)[names(col.vec)%in%host.legend$old.colors],
                    labels=host.legend$new.colors)+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", paste(lvl.name,collapse = " and ")))+
  theme(plot.margin=unit(c(1,1,1,2.5), 'cm'),
        strip.text.x = ggtext::element_markdown(size = 30),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=30,hjust=1),
        axis.text.y = element_text(size=30),
        axis.title = element_text(size = 30),
        plot.title = ggtext::element_markdown(size = 30),
        plot.caption = element_text(size=23),
        legend.text = element_markdown(size = 30),
        legend.title = element_text(size = 35),
        legend.position = "right")


for(image.format in image.formats){
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "barplot","NMR-B6mouse",
                      truncationlvl,agglom.rank,
                      sep = "-"),".",image.format),
         plot=lvl.plot,
         width = 9000,height = 6000,
         units = "px",dpi=300,device = image.format)
}
## Barplot for NMR and NMR wt ####
if("NMRwt"%in%custom.levels){
  lvl.df<-ps.q.agg.for_bp%>%
    group_by(Sample)%>%
    mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
    # so we can have info about sample size
    ungroup()%>%
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
    filter(class%in%c("NMR","NMRwt"),Abundance!=0)
  lvl.name<-unname(pretty.level.names[names(pretty.level.names)%in%c("NMR","NMRwt")])
  lvl.name<-gsub("<br>"," ", lvl.name)
  # the total legend is big, we need to narrow down to our host
  host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon.bp%in%names(table(lvl.df$Taxon.bp))]
  nmr.set<-lvl.df%>%
    filter(class=="NMR")%>%
    select(Taxon.bp)%>%
    unique()%>%
    pull()
  
  nmrwt.set<-lvl.df%>%
    filter(class=="NMRwt")%>%
    select(Taxon.bp)%>%
    unique()%>%
    pull()
  
  nmr.uniq<-setdiff(nmr.set,nmrwt.set)
  nmr.uniq.legend<-ps.q.legend$Taxon.bp[ps.q.legend$Taxon.bp%in%nmr.uniq]
  
  # New font colors
  host.legend<-data.frame(host.legend,host.legend)%>%
    rename(old.colors="host.legend",
           new.colors="host.legend.1")%>%
    mutate(new.colors=ifelse(host.legend%in%nmr.uniq.legend,
                             paste("<span style='color: red'>",new.colors,"</span>"),
                             old.colors))
  
  lvl.plot<-lvl.df%>%
    ggplot(aes(x=NewSample, y=RelativeAbundance,  
               fill=factor(Taxon.bp, levels=host.legend$old.colors)))+
    # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
    geom_bar(stat = "identity")+ # barplot
    facet_grid(~class, # separate species
               scales="free",  # each species will have its own bars inside facet (instead of all bars)
               space = "free", # bars will have same widths
               labeller = labeller(class=pretty.level.names) # labeller will change facet labels to custom
    )+
    guides(fill=guide_legend(ncol=1))+ # legend as one column
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values = col.vec[names(col.vec)%in%host.legend$old.colors],
                      breaks = names(col.vec)[names(col.vec)%in%host.legend$old.colors],
                      labels=host.legend$new.colors)+
    xlab("") +
    ylab("Relative Abundance (%)")+
    labs(fill="Taxon")+
    theme_bw()+
    ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", paste(lvl.name,collapse = " and ")))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          strip.text.x = ggtext::element_markdown(size = 20),
          panel.spacing = unit(0.8, "cm"), # increase distance between facets
          axis.text.x = element_text(angle=45,size=20,hjust=1),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = ggtext::element_markdown(size = 25),
          plot.caption = element_text(size=23),
          legend.text = element_markdown(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  for(image.format in image.formats){
    ggsave(paste0(barplot.directory,
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "barplot","NMR-NMRwt",
                        truncationlvl,agglom.rank,
                        sep = "-"),".",image.format),
           plot=lvl.plot,
           width = 9000,height = 6000,
           units = "px",dpi=300,device = image.format)
  }
}

# Session Info
sessionInfo()


# Average barplots ####
# For NMR
ps.q.agg%>%
  filter(class=="NMR",MeanRelativeAbundance>=1)%>%
  group_by_at(agglom.rank)%>%
  distinct(.,get(agglom.rank),.keep_all = TRUE)%>%
  ungroup()%>%
  dplyr::select(Phylum,all_of(agglom.rank),class,MeanRelativeAbundance)%>%
  # pivot_wider(names_from = class,
              # values_from = MeanRelativeAbundance,
              # values_fill = 0)%>%
  dplyr::arrange(-MeanRelativeAbundance)%>%
  dplyr::mutate(csumNMR=cumsum(MeanRelativeAbundance))

ps.q.agg%>%
  filter(class=="B6mouse",MeanRelativeAbundance>=1)%>%
  group_by_at(agglom.rank)%>%
  distinct(.,get(agglom.rank),.keep_all = TRUE)%>%
  ungroup()%>%
  dplyr::select(Phylum,all_of(agglom.rank),class,MeanRelativeAbundance)%>%
  # pivot_wider(names_from = class,
  # values_from = MeanRelativeAbundance,
  # values_fill = 0)%>%
  dplyr::arrange(-MeanRelativeAbundance)%>%
  dplyr::mutate(csumNMR=cumsum(MeanRelativeAbundance))

# barplot of agglom.rank taxa (e.g. genera) with average abundance >=1% in NMR
ps.q.agg%>%
  filter(class=="NMR",MeanRelativeAbundance>=1)%>%
  distinct(get(agglom.rank),.keep_all = TRUE)%>%
  ggplot(.,aes(x=get(agglom.rank),y=MeanRelativeAbundance))+
  geom_bar(stat = "identity")+
  coord_flip()


####
# Create a barplot of relative abundances for a certain taxon across hosts
# (quick and dirty)
ggplot.taxon<-function(tax.df,taxon.to.plot,tax.rank){
  ggplot.obj<-tax.df%>%
    filter(get(tax.rank)==taxon.to.plot,
           Abundance!=0)%>%
    ggplot(aes(x=Sample,y=RelativeAbundance),fill=tax.rank)+
    geom_bar(stat="identity")+
    facet_grid(.~class,scales = "free",
               space = "free")+
    theme_bw()+
    labs(y="Relative abundance (%)",
         title = paste(taxon.to.plot, "relative abundance across rodents"))+
    theme(axis.text.x = element_text(angle=45,hjust=1))
  return(ggplot.obj)
}

ggplot.taxon(ps.q.agg,"Lactobacillus","Genus")

# ggsave(filename = paste0("./images/barplots/Lactobacillus.png"),
#        last_plot(),
#        height = 800,width = 1500,units = "px",dpi = 120)
ggplot.taxon(ps.q.agg,"Bifidobacterium","Genus")

# ggsave(filename = paste0("./images/barplots/Bifidobacterium.png"),
#        last_plot(),
#        height = 800,width = 1500,units = "px",dpi = 120)
ggplot.taxon(ps.q.agg,"Desulfovibrio","Genus")
# ggsave(filename = paste0("./images/barplots/Desulfovibrio.png"),
#        last_plot(),
#        height = 800,width = 1500,units = "px",dpi = 120)
