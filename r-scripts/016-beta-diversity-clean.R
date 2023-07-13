library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)
library(pals)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"
# load("./rdafiles/yashliu-dada2-284-203-Genus-138-1-phyloseq-workspace.RData")

load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

pretty.facet.labels<-
  c("NMR" = "*Heterocephalus glaber*", # better labels for facets
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
custom.levels<-names(pretty.facet.labels)

# Using Polychrome package
set.seed(1)
# custom.colors<- createPalette(length(custom.levels),
#               seedcolors = rainbow(6))
custom.colors<- 
  createPalette(length(custom.levels),
                seedcolors = c("#B22222", "#0000FF","#006400", 
                               "#FF8C00", "#68228B"))
custom.colors<-unname(custom.colors)
swatch(custom.colors)
# custom colors for scale_fill_manual (maybe not needed)
# custom.fill<-c("dodgerblue", "seagreen","indianred")

# custom labels for scale_color_manual
custom.color.labels<-unname(pretty.facet.labels)

# Choose distance metric
dist.metric<-"bray"
dist.metric<-"jaccard"
dist.metric<-"aitchison"
dist.metric<-"canberra"
if(dist.metric=="bray"){
  beta.label<-"Bray-Curtis dissimilarities"
}else {
  beta.label<-paste(str_to_title(dist.metric),"distances")
}

permut.num<-1000 # number of permutations for PERMANOVA

ps.sampledata<-ps.q.agg.abs%>%
  select(c("Sample","class","origin","relation",
           "sex","caste","birthday"))%>%
  distinct() # metadata

ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)

ps.q.df <-ps.q.agg.rel%>%
  select(Sample,OTU,Abundance,class,Taxon)
# ps.q.df <-ps.q.agg.abs%>%
#   select(Sample,OTU,Abundance,class)

# Beta diversity ####
# find the smallest sample size
min.n_seqs.all<-ps.q.df%>%
  select(Sample,OTU,Abundance)%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# rarefaction
# ps.rarefied = rarefy_even_depth(ps.q.df,
#                                 sample.size=17000, 
#                                 replace=F)

# convert the data frame into wide format
ps.q.df.wide<-ps.q.df%>%
  select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  
ps.q.df.wide<-as.matrix(ps.q.df.wide)


# Bray or Jaccard ####
# Bray is meaningful only for integers (counts)
set.seed(1)
# Rarefaction is done by avgdist
# computes the dissimilarity matrix of a dataset MULTIPLE times using vegdist
# while randomly subsampling the dataset each time.
# All of the subsampled iterations are then averaged (mean) to provide a 
# distance matrix that represents the average of multiple subsampling iterations.

if(dist.metric=="jaccard"){
  dist<-avgdist(ps.q.df.wide,
                dmethod="jaccard",
                binary=TRUE,
                sample=min.n_seqs.all,
                iterations = 1000)
}else if(dist.metric=="bray"|dist.metric=="canberra"){
  dist<-avgdist(ps.q.df.wide,
                dmethod=dist.metric,
                sample=min.n_seqs.all,
                iterations = 1000)
}else if(dist.metric=="aitchison"){
  # dist.tfm<-decostand(x = ps.q.df.wide,"clr",pseudocount=1)
  dist<-avgdist(ps.q.df.wide,
                dmethod="aitchison",
                sample=min.n_seqs.all,
                iterations = 1000,
                pseudocount=1)
}
set.seed(1)


dist.df<-dist%>% # convert distances into tibble
  as.matrix()%>%
  as_tibble()
dist.df<-cbind(labels(dist),dist.df) # create a column of sample ids
colnames(dist.df)[1]<-"Sample"

meta.dist<-dist.df%>%inner_join(ps.sampledata,.,by="Sample") # distances with metadata

all.dist<-meta.dist%>%
  select(all_of(.$Sample))%>% # extract distances (corresponding to Sample vector)
  as.dist()

## PERMANOVA ####
all.test<-adonis2(all.dist~class, # distances explained by class
                  data=meta.dist, 
                  permutations = permut.num) # permutation number
all.test$F[1]
all.test.p<-all.test$`Pr(>F)`[1]
all.test.p

## Making pairwise comparisons ####
adonis2(all.dist~class,
        data=meta.dist,
        permutations = permut.num)
str(meta.dist)

# we found a significant p value, but which groups are different?
pairwise_p<-numeric()
pairwise_F<-numeric()

# pairwise tests
combinations<-combn(custom.levels,2)
for (i in 1:ncol(combinations)) {
  test.data<-meta.dist%>%
    filter(class==combinations[1,i] | class==combinations[2,i]) # extract data for two groups
  test.distances<-test.data%>%# extract only distances
    select(all_of(.$Sample))%>% # extract distances (corresponding to Sample vector)
    as.dist()
  adonis.test<-adonis2(test.distances~class,
                       data=test.data,
                       permutations=permut.num)
  pairwise_p[paste(combinations[1,i],combinations[2,i],sep = "_")]<-
    adonis.test$`Pr(>F)`[1]
  pairwise_F[paste(combinations[1,i],combinations[2,i],sep = "_")]<-
    adonis.test$F[1]
}

### Results ####
pairwise_p
pairwise_F

### adjusting p values ####
p.adjust(pairwise_p, method = "BH")
which(pairwise_p>0.05)
stopifnot(all(p.adjust(pairwise_p, method = "BH")<0.05))

## NMDS ####
# performs Nonmetric Multidimensional Scaling (NMDS), and tries to find 
# a stable solution using several random starts. 
# In addition, it standardizes the scaling in the result, 
# so that the configurations are easier to interpret, and adds species scores 
# to the site ordination. The metaMDS function does not provide actual NMDS, 
# but it calls another function for the purpose. 
# Currently monoMDS is the default choice, and it is also possible to call
# the isoMDS (MASS package).
set.seed(1)
nmds<-metaMDS(dist, autotransform = FALSE) 
# autotransform = FALSE means NO Square root transformation
# and NO Wisconsin double standardization

nmds.scores<-scores(nmds)%>%
  as_tibble(rownames="Sample")%>% # convert rownames into "Sample" column
  inner_join(.,ps.sampledata, by="Sample") 


nmds.plot<-ggplot(nmds.scores,
                    aes(x=NMDS1,y=NMDS2,color=class,fill=class))+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  # scale_color_manual(breaks = unique(sort(sample.lookup$host)),
  #                    labels=unique(sort(sample.lookup$host)),
  theme_bw()+
  guides(fill="none")+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  labs(color="Host"
       # ,
       # caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 using Benjamini-Hochberg correction for multiple comparisons"
  )+
  # ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE)+ # add labels to samples
  ggtitle(paste0("nMDS on ", beta.label, " between different rodents (",agglom.rank, " level)"))+
  theme(
    axis.text.x = element_text(angle=0,size=20,hjust=1),
    axis.text.y = element_text(size=20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 27),
    legend.text = ggtext::element_markdown(size = 20),
    legend.title = element_text(size = 25),
    legend.position = "right",
    plot.caption = ggtext::element_markdown(hjust = 0, size=20),
    plot.caption.position = "plot")

ggsave(paste0("./images/diversity/",Sys.Date(),"-ndms-",dist.metric,"-all-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=nmds.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/",Sys.Date(),"-nmds-",dist.metric,"-all-",agglom.rank,
              "-",truncationlvl,".tiff"),plot=nmds.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "tiff")


# rarefy without replacement to the min sample size
# set.seed(1)
# ps.q.agg.abs.rarefied = rarefy_even_depth(ps.q.agg.abs, 
#                                         rngseed=1, 
#                                         sample.size=min(sample_sums(ps.q.agg.abs)), 
#                                         replace=F)
# ps
# ps.q.agg.abs.rarefied
# plot_bar(ps.q.agg.abs.rarefied, fill="Genus")+
#   facet_grid(~animal, scales="free",space="free")+  
#   guides(fill=guide_legend(ncol=1))

## PCoA ####
## Bray or Jaccard
pcoa<-cmdscale(dist,
                 k=2,
                 eig = TRUE,
                 add = TRUE) # PCoA
# k is the num of principal coordinates
# eig allows to see % of variation explained
# add rescales eigenvalues to make them all positive
positions<-pcoa$points # pcoa values to plot
colnames(positions)<-c("pcoa1", "pcoa2")

percent_explained<-round(100* pcoa$eig / 
                             sum(pcoa$eig),1)
# previous won't have 0 after decimal

# format() command will show exact number of digits you need 
percent_exp<-format(round(100* pcoa$eig / 
                              sum(pcoa$eig),1),1)

# % explained by each axis is the value in eig divided 
# by the sum of eig and multiplied by 100
# we create a vector of % (each axis is divided by sum(eig))
# so, we actually have % for all axes
# first two axes are percent_explained[1:2]
pcoa.positions<-positions %>% 
  as_tibble(rownames="Sample") %>%
  inner_join(.,ps.sampledata, by="Sample")  

pcoa.plot<-ggplot(pcoa.positions,
                    aes(x=pcoa1,y=pcoa2,color=class,
                        fill=class))+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  theme_bw()+
  guides(fill="none")+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  labs(x=paste0("PCo 1 (", percent_exp[1],"%)"),
       y=paste0("PCo 2 (", percent_exp[2],"%)"),
       color="Host"
       # ,
       # caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 using Benjamini-Hochberg correction for multiple comparisons"
  )+
  # ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE)+ # add labels to samples
  ggtitle(paste0("PCoA on ", beta.label, " between different rodents (",agglom.rank, " level)"))+
  theme(
    axis.text.x = element_text(angle=0,size=20,hjust=1),
    axis.text.y = element_text(size=20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 27),
    legend.text = ggtext::element_markdown(size = 20),
    legend.title = element_text(size = 25),
    legend.position = "right",
    plot.caption = ggtext::element_markdown(hjust = 0, size=20),
    plot.caption.position = "plot")



ggsave(paste0("./images/diversity/",Sys.Date(),"-pcoa-",dist.metric,"-all-",
              agglom.rank, "-",truncationlvl,".png"),
       plot=pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")

ggsave(paste0("./images/diversity/",Sys.Date(),"-pcoa-",dist.metric,"-all-",
              agglom.rank, "-",truncationlvl,".tiff"),
       plot=pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "tiff")
rm(list=ls())


# Unweighted UniFrac ####
# Use relative abundances
ps.q.df <-ps.q.agg.rel%>%
  select(Sample,OTU,Abundance,class,Taxon)
ps.q.df.wide<-ps.q.df%>%
  select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  
ps.q.df.wide<-as.matrix(ps.q.df.wide)

uni.unw<-UniFrac(ps.q, weighted = FALSE)
uni.unw.df<-uni.unw%>% # convert distances into tibble
  as.matrix()%>%
  as_tibble()

uni.unw.df<-cbind(labels(uni.unw),uni.unw.df) # create a column of sample ids
colnames(uni.unw.df)[1]<-"Sample"

## Unweighted unifrac PCoA ####
uni.unw.pcoa<-cmdscale(uni.unw,
                       k=2,
                       eig = TRUE,
                       add = TRUE)

uni.unw.positions<-uni.unw.pcoa$points # pcoa values to plot
colnames(uni.unw.positions)<-c("pcoa1", "pcoa2")

uni.unw.percent_explained<-round(100* uni.unw.pcoa$eig / 
                                   sum(uni.unw.pcoa$eig),1)
# previous won't have 0 after decimal

# format() command will show exact number of digits you need 
uni.unw.percent_exp<-format(round(100* uni.unw.pcoa$eig / 
                                    sum(uni.unw.pcoa$eig),1),1)

# % explained by each axis is the value in eig divided 
# by the sum of eig and multiplied by 100
# we create a vector of % (each axis is divided by sum(eig))
# so, we actually have % for all axes
# first two axes are percent_explained[1:2]
uni.unw.pcoa.positions<-uni.unw.positions %>% 
  as_tibble(rownames="Sample") %>%
  inner_join(.,ps.sampledata, by="Sample")  

ggplot(uni.unw.pcoa.positions,
       aes(x=pcoa1,y=pcoa2,color=class,
           fill=class))+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  theme_bw()+
  guides(fill="none")+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)


uni.unw.pcoa.plot<-ggplot(data = uni.unw.pcoa.positions, 
                          aes(x = pcoa1, y = pcoa2,color = class)) + 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  # geom_point(size = 3, alpha = 0.5)+
  geom_point()+
  theme_bw() + 
  guides(fill="none")+
  scale_colour_manual(breaks = custom.levels,
                      labels=custom.color.labels,
                      values = custom.colors) + 
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  ggtitle(paste0("PCoA on unweighted UniFrac distances between different rodents
                 (",agglom.rank, " level)"))+
  labs(x=paste0("PCo 1 (", uni.unw.percent_exp[1],"%)"),
       y=paste0("PCo 2 (", uni.unw.percent_exp[2],"%)"),
       color = "Host")+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 25), 
        legend.text = ggtext::element_markdown(size = 20),
        legend.position = "right",
        plot.caption = ggtext::element_markdown(hjust = 0, size=20),
        plot.caption.position = "plot")


ggsave(paste0("./images/diversity/",Sys.Date(),"-uni-unw-pcoa-all-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=uni.unw.pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/",Sys.Date(),"-uni-unw-pcoa-all-",agglom.rank,
              "-",truncationlvl,
              ".tiff"),plot=uni.unw.pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "tiff")
