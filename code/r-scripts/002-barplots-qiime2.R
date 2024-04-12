# Creating barplots ####

# After importing the QZA files into R and processing them with phyloseq,
# it's time to explore the taxonomic composition of our data.
# We will use the Polychrome package to create a custom palette for the 
# barplots.
library(phyloseq)
library(tidyverse)
library(Polychrome)
library(ggtext)
## Specifying parameters and directory/file names #### 
authorname<-"pooled" # name of the folder with QIIME2 output
agglom.rank<-"Order" # this is the taxonomic rank that was used for agglomeration
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


# Ordering the legend
if(asvlevel==TRUE){
  taxa.list<-ps.q.agg%>%
    group_by(class,OTU)%>%
    filter(MeanRelativeAbundance>=0.1)%>%
    ungroup()%>%
    select(matches(paste0("^",agglom.rank,"$")))%>%
    pull(.)%>%
    unique()
}else{
  taxa.list<-ps.q.agg%>%
    group_by_at(c("class",agglom.rank))%>%
    filter(MeanRelativeAbundance>=0.1)%>%
    ungroup()%>%
    select(matches(paste0("^",agglom.rank,"$")))%>%
    pull(.)%>%
    unique()
}

taxa.list<-c(taxa.list,"Remainder")
custom_order <- c("Remainder","Kingdom", "Phylum", "Class", "Order", "Family")
# If we agglomerate by higher level (Order,Class, etc), need to adjust the rank
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}

# Get only classified taxa to create Genus + Family names
rank <- match(sub(".* ", "", taxa.list),custom_order)
rank[is.na(rank)] <- length(custom_order) + 1

# Only agglom.rank: their rank is the last in the custom_order vector
classified.taxa<-taxa.list[rank==length(custom_order) + 1]
unclassified.taxa<-taxa.list[rank!=length(custom_order) + 1]

if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}

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
# Order by higher rank then agglom.rank
# These are agglom.rank (preceding.rank) format strings for barplot
taxa.for_bp.df<-taxa.for_bp.df%>%
  arrange(get(preceding.rank),get(agglom.rank))
taxa.for_bp.list<-taxa.for_bp.df$Taxon.bp

# Now we order unclassified taxa
newrank <- match(sub(".* ", "", unclassified.taxa),custom_order)
# newrank[is.na(newrank)] <- length(custom_order) + 1

unclassified.taxa<-unclassified.taxa[order(newrank)]

unclassified.taxa.split <- split(unclassified.taxa, newrank[order(newrank)])
unclassified.taxa.sorted <- unlist(lapply(unclassified.taxa.split, sort))
unclassified.taxa.sorted<-unname(unclassified.taxa.sorted)

# Add sorted unclassified taxa to the sorted classified taxa
taxa.for_bp.list<-c(unclassified.taxa.sorted,taxa.for_bp.list)
taxa.for_bp.list[1]<-"Remainder (Mean abundance < 0.1%)"
# Add barplot taxa to the main dataframe
ps.q.agg<-ps.q.agg%>%
  left_join(taxa.for_bp.df,by=agglom.rank)%>%
  ungroup()%>%
  select(-paste0(preceding.rank,".y"))%>% # remove the preceding.rank column from taxa.for_bp.df
  rename(!!preceding.rank:=paste0(preceding.rank,".x")) # rename the preceding.rank column
# !!preceding.rank:= will evaluate the variable

ps.q.agg[which(ps.q.agg$MeanRelativeAbundance>=0.1&is.na(ps.q.agg$Taxon.bp)),"Taxon.bp"]<-
  ps.q.agg[which(ps.q.agg$MeanRelativeAbundance>=0.1&is.na(ps.q.agg$Taxon.bp)),agglom.rank]

ps.q.agg[which(ps.q.agg$MeanRelativeAbundance<0.1),"Taxon.bp"]<-
  "Remainder (Mean abundance < 0.1%)"

# We want to highlight NMR-specific taxa
nmr.set<-ps.q.agg%>%
  filter(class=="NMR")%>%
  select(all_of(agglom.rank))%>%
  unique()%>%
  pull()

others.set<-ps.q.agg%>%
  filter(class!="NMR")%>%
  select(all_of(agglom.rank))%>%
  unique()%>%
  pull()

nmr.uniq<-setdiff(nmr.set,others.set)
agglom.rank.vec<-ps.q.agg%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=0.1)%>%
  select(agglom.rank)%>%
  pull()%>%
  unique()
nmr.uniq.legend<-agglom.rank.vec[agglom.rank.vec%in%nmr.uniq]

# New font colors
ps.q.legend<-as.data.frame(taxa.for_bp.list)%>%
  rename("Taxon.bp"="taxa.for_bp.list")%>%
  mutate(new.colors=ifelse(Taxon.bp%in%nmr.uniq.legend,
                           paste("<span style='color: red'><b>",Taxon.bp,"</b></span>"),
                           Taxon.bp))
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
plot.cols<-createPalette(length(taxa.for_bp.list)-1,
                         seedcolors =rainbow(7))# input: number of rows

# in our legend and the seed colors that we decide to be rainbow

# The vector of colors should be named according to our legend
# because we will use this mapping in the barplot. So, each color in the vector
# correspond to each color in the legend. Remember, remainder taxa are 
# all merged into a single entry ("Remainder"), so there's just one color for 
# remainder portion.

# Decrease alpha for unclassified
plot.cols<-c("#C1CDCD",plot.cols)
plot.cols[which(taxa.for_bp.list%in%unclassified.taxa.sorted)]<-
  adjustcolor(plot.cols[which(taxa.for_bp.list%in%unclassified.taxa.sorted)],
              alpha.f = 0.5)

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
mainplot<-ps.q.agg%>%
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
             labeller = labeller(class=pretty.facet.labels))+ # labeller will
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
        legend.text = element_markdown(size = 20), # size of legend text
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
  lvl.name<-unname(pretty.facet.labels[i]) # We find the pretty name for the 
  # facet using the `pretty.facet.labels` vector. `unname` will remove
  # the name from the vector element (name was taken from custom.levels, 
  # not pretty)
  lvl.name<-gsub("<br>"," ", lvl.name) # also remove all line breaks
  # the total legend is big, we need to narrow down to our host. 
  # Take the legend and extract taxa that are present in the lvl.df
  host.legend<-ps.q.legend$Taxon.bp[ps.q.legend$Taxon.bp%in%names(table(lvl.df$Taxon.bp))]
  
  
  lvl.plot<-lvl.df%>%
    ggplot(aes(x=NewSample, y=RelativeAbundance,  
               fill=factor(Taxon.bp, levels=host.legend)))+
    # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
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
  group_by(Sample)%>%
  mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
  # so we can have info about sample size
  ungroup()%>%
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of 
  # our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
  filter(class%in%c("NMR","B6mouse"),Abundance!=0)
lvl.name<-unname(pretty.facet.labels[names(pretty.facet.labels)%in%c("NMR","B6mouse")])
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
nmr.uniq.legend<-ps.q.agg%>%
  filter(Taxon.bp%in%taxa.for_bp.list,MeanRelativeAbundance>=0.1)%>%
  distinct(get(agglom.rank),Taxon.bp)%>%
  rename(!!agglom.rank:="get(agglom.rank)")%>%
  filter(get(agglom.rank)%in%nmr.uniq)%>%
  select(Taxon.bp)%>%
  pull()

# New font colors
host.legend<-data.frame(host.legend,host.legend)%>%
  rename(old.colors="host.legend",
         new.colors="host.legend.1")%>%
  mutate(new.colors=ifelse(host.legend%in%nmr.uniq.legend,
                           paste("<span style='color: red'><b>",new.colors,"</b></span>"),
                           old.colors))

lvl.plot<-lvl.df%>%
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=host.legend$old.colors)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
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
        axis.line = element_blank(), 
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
if("NMRwt"%in%custom.levels){
  lvl.df<-ps.q.agg%>%
    group_by(Sample)%>%
    mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
    # so we can have info about sample size
    ungroup()%>%
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
    filter(class%in%c("NMR","NMRwt"),Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[names(pretty.facet.labels)%in%c("NMR","NMRwt")])
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
               labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
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
          axis.line = element_blank(), 
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
}

# Session Info
sessionInfo()

plot1<-ps.q.agg%>%
  filter(Genus=="Lactobacillus",
         Abundance!=0)%>%
  ggplot(aes(x=Sample,y=RelativeAbundance),fill=Genus)+
  geom_bar(stat="identity")+
  facet_grid(.~class,scales = "free",
             space = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Lactobacillus")
ggsave(filename = paste0("./images/barplots/Lactobacillus.png"),
       plot1,
       height = 800,width = 1500,units = "px",dpi = 120)

plot2<-ps.q.agg%>%
  filter(Genus=="Bifidobacterium",
         Abundance!=0)%>%
  ggplot(aes(x=Sample,y=RelativeAbundance),fill=Genus)+
  geom_bar(stat="identity")+
  facet_grid(~class,scales = "free",
             space = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Bifidobacterium")
ggsave(filename = paste0("./images/barplots/Bifidobacterium.png"),
       plot2,
       height = 800,width = 1500,units = "px",dpi = 120)
plot3<-ps.q.agg%>%
  filter(Genus=="Desulfovibrio",
         Abundance!=0)%>%
  ggplot(aes(x=Sample,y=RelativeAbundance),fill=Genus)+
  geom_bar(stat="identity")+
  facet_grid(~class,scales = "free",
             space = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Desulfovibrio")
ggsave(filename = paste0("./images/barplots/Desulfovibrio.png"),
       plot3,
       height = 800,width = 1500,units = "px",dpi = 120)
