authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"

# Import data ####
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
# ps.q.df.preprocessed.date_time is a rarefied dataset for all hosts
ps.q.df.preprocessed.date_time<-"20240426_22_00_04" # all hosts, genus level
image.formats<-c("png","tiff")

pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "NMRwt"="Wild *Heterocephalus glaber*",
                      "DMR" = "*Fukomys Damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*")

excluded.samples<-
  c(#"MSMmouse",
    #"FVBNmouse",
    "NMRwt")
custom.levels<-names(pretty.level.names)

# filter your data
custom.levels<-custom.levels[!custom.levels%in%excluded.samples]
# load the output of 003-phyloseq-rarefaction-filtering.R file ####
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)
##########################

barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved
image.formats<-c("png","tiff")
# phyloseq.workspace.date_time is an Rdata workspace with nonrarefied dataframe
# from 001-phyloseq-qiime2.R
phyloseq.workspace.date_time="20240426_21_44_30" 
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

#######################
custom.md.ages<-readRDS("./output/rdafiles/custom.md.ages.rds")
custom.md.ages<-custom.md.ages%>%
  mutate(old_agegroup=gsub("\\[0\\,10\\)","Young group",old_agegroup),
         old_agegroup=gsub("\\[10\\,16\\)","Old group",old_agegroup))
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  select(Sample,agegroup,old_agegroup)%>%
  arrange(agegroup)%>%
  distinct()
# pretty.agegroup.names<-names(table(custom.md.ages$old_agegroup))
# names(pretty.agegroup.names)<-names(table(custom.md.ages$agegroup))
pretty.agegroup.names<-c("agegroup0_10"="Young group",
                         "agegroup10_16"="Old group")
agegroup.levels<-names(pretty.agegroup.names)


# Non rarefied df plot
lvl.df%>%
  left_join(custom.md.ages)%>%
  # mutate(Sample=factor(Sample,levels=sample.levels))
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=host.legend)))+
  geom_bar(stat = "identity")+ # barplot
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  facet_grid(~agegroup, # separate species
             scales="free_x",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(agegroup=pretty.agegroup.names) )+# labeller will change facet labels to custom
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


# 
# ggsave(paste0(barplot.directory,
#               paste(paste(format(Sys.time(),format="%Y%m%d"),
#                           format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                     custom.levels[i],"barplot",truncationlvl,
#                     agglom.rank,sep = "-"),".png"),
#        plot=lvl.plot,
#        width = 8000,height = 6000,
#        units = "px",dpi=300,device = "png")
colnames(lvl.df)
colnames(ps.q.df.preprocessed)




ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
  filter(class=="NMR")%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample",agglom.rank))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)

ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class",agglom.rank))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(-TotalClass,-TotalAgglomRank)

ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
  left_join(taxa.for_bp.df)

ps.q.df.preprocessed[which(ps.q.df.preprocessed$MeanRelativeAbundance>=1&is.na(ps.q.df.preprocessed$Taxon.bp)),"Taxon.bp"]<-
  ps.q.df.preprocessed[which(ps.q.df.preprocessed$MeanRelativeAbundance>=1&is.na(ps.q.df.preprocessed$Taxon.bp)),agglom.rank]
# Taxa with MeanRelativeAbundance<1% become "Remainder"
ps.q.df.preprocessed[which(ps.q.df.preprocessed$MeanRelativeAbundance<1),"Taxon.bp"]<-
  "Remainder (Mean relative abundance < 1%)"

ps.rare<-ps.q.df.preprocessed%>%
  ungroup()%>%
  select(Genus,RelativeAbundance,Sample)%>%
  arrange(Sample,-RelativeAbundance)
ps.nonrare<-lvl.df%>%
  ungroup()%>%
  select(Genus,RelativeAbundance,Sample)%>%
  arrange(Sample,-RelativeAbundance)

table(ps.rare$Genus%in%ps.nonrare$Genus)
ps.nonrare[!ps.nonrare$Genus%in%ps.rare$Genus,]

colnames(ps.nonrare)<-c("Genus","RelativeAbundance_nr","Sample")
ps.rare%>%
  full_join(ps.nonrare,by=c("Sample"="Sample",
                         "Genus"="Genus"))%>%
  mutate(diff.ra=abs(RelativeAbundance-RelativeAbundance_nr))%>%
  arrange(-diff.ra)%>%View


ps.q.df.preprocessed%>%
  select(-TotalSample)%>%
  group_by(Sample)%>%
  select(-birthday)%>%
  mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
  # so we can have info about sample size
  ungroup()%>%
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of
  # our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ", TotalAbundance, ")"))%>% # add a 
  # column where sample names are together with sample sizes
  left_join(custom.md.ages)%>%
  # mutate(Sample=factor(Sample,levels=sample.levels))
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=host.legend)))+
  geom_bar(stat = "identity")+ # barplot
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  facet_grid(~agegroup, # separate species
             scales="free_x",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(agegroup=pretty.agegroup.names) )+# labeller will change facet labels to custom
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
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
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
                    custom.levels[i],"barplot-rare",truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=last_plot(),
       width = 8000,height = 6000,
       units = "px",dpi=300,device = "png")
