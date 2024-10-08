# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan","Polychrome"))
library(tidyverse)
library(vegan)
library(Polychrome)
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"

# Import data ####
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
# Import abundance table from 001-phyloseq-qiime2.R
ps.q.agg.date_time<-"20240620_12_40_41"

ps.q.agg<-readRDS(file=file.path(
  "./output/rdafiles",
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))

# 20240620_12_38_18 ps.q.agg.date_time OTU level, all hosts
# 20240620_12_40_41 ps.q.agg.date_time genus level, all hosts
# 20240917_10_54_12 ps.q.agg.date_time family level, all hosts
# 20240917_21_29_36 ps.q.agg.date_time phylum level, all hosts

ps.q.df.pca.input.date_time<-"20240426_22_00_04"
# Import metadata
custom.md<-readRDS("./output/rdafiles/custom.md.rds")

image.formats<-c("png","tiff")
plot.types<-c("plot"="",
              "plot.with.labels"="-with-labels")

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
  c("NMRwt")

custom.levels<-intersect(names(pretty.level.names),custom.md$class)

# Create custom palette with Polychrome package
set.seed(1)
custom.colors<- 
  createPalette(length(custom.levels),
                seedcolors = c("#B22222","#5D478B", "#0000FF","#006400", 
                               "#FF8C00","#00EE00","#EEC900", "#00FFFF","#FF6EB4"))
# custom.colors<-unname(custom.colors)
names(custom.colors)<-custom.levels
swatch(custom.colors)

# Setup general theme for ggplot2. Add more parameters depending on the plot
mytheme <- theme(plot.title = element_text(size = 27),
                 axis.text.x = element_text(angle=0,size=20,hjust=1),
                 axis.text.y = element_text(size=20),
                 axis.title = element_text(size = 20),
                 legend.text = ggtext::element_markdown(size = 20),
                 legend.title = element_text(size = 25),
                 legend.position = "right",
                 plot.caption = ggtext::element_markdown(hjust = 0, size=20),
                 plot.caption.position = "plot"
)
# filter your data
if(exists("excluded.samples")){
  custom.levels<-custom.levels[!custom.levels%in%excluded.samples]
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  pretty.level.names<-pretty.level.names[!names(pretty.level.names)%in%excluded.samples]
  ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,!class%in%excluded.samples,Abundance!=0)
}else {
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
}

if(exists("excluded.samples")){
  custom.colors<-custom.colors[!names(custom.colors)%in%excluded.samples]
}
# Choose distance metric
dist.metric<-"robust.aitchison"
dist.metric<-"jaccard"
dist.metric<-"canberra"
dist.metric<-"bray"
if(dist.metric=="bray"){
  beta.label<-"Bray-Curtis dissimilarities"
}else { # in case of robust.aitchison, we need to substitute the dot
  beta.label<-paste(str_to_title(gsub("\\.", " ", dist.metric)),"distances")
}

permut.num<-1000 # number of permutations for PERMANOVA

# Beta diversity ####
# find the smallest sample size
min.n_seqs.all<-ps.q.agg%>%
  select(Sample,Abundance,class,all_of(agglom.rank))%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# convert the data frame into wide format
ps.q.df.wide<-ps.q.agg%>%
  select(Sample,Abundance,class,all_of(agglom.rank))%>%
  pivot_wider(names_from = all_of(agglom.rank), # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs/genera/species and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
# Remove the first two columns (Sample and class)
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]
# Convert the dataframe into matrix
ps.q.df.wide<-as.matrix(ps.q.df.wide)


# Calculate distances ####
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
}else if(dist.metric=="robust.aitchison"){
  dist<-avgdist(ps.q.df.wide,
                dmethod="robust.aitchison",
                sample=min.n_seqs.all,
                iterations = 1000
                )
}

dist.df<-dist%>% # convert distances into tibble
  as.matrix()%>%
  as_tibble()
dist.df<-cbind(labels(dist),dist.df) # create a column of sample ids
colnames(dist.df)[1]<-"Sample"

meta.dist<-dist.df%>%inner_join(custom.md,.,by="Sample") # distances with metadata

## PERMANOVA ####
# To test whether there were significant differences in beta diversity between 
# animal hosts, we perform PERMANOVA with the adonis2 function
all.test<-adonis2(dist~class, # distances explained by class
                  data=meta.dist, 
                  permutations = permut.num) # permutation number
all.test$F[1]
all.test.p<-all.test$`Pr(>F)`[1]
all.test.p

## Making pairwise comparisons ####
adonis2(dist~class,
        data=meta.dist,
        permutations = permut.num)
str(meta.dist)

# we found a significant p value, but which groups are different?
pairwise_p<-numeric()
pairwise_F<-numeric()

# pairwise tests
combinations<-combn(custom.levels,2)
for (i in 1:ncol(combinations)) {
  test.data<-meta.dist%>% # combinations is a matrix 2xn matrix, so 
    # each combination is stored in a column (two rows, so it's a pair)
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
stopifnot(all(p.adjust(pairwise_p, method = "BH")<0.05)) # B6mouse_FVBNmouse 

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
  inner_join(.,custom.md, by="Sample") 


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
                     labels=unname(pretty.level.names),
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  labs(color="Host")+
  ggtitle(paste0("nMDS on ", beta.label, " between different rodents (",agglom.rank, " level)"))+
  mytheme # general theme


# Depending on the results of PERMANOVA, we should choose an appropriate
# label
if(all.test.p<0.05 &sum(pairwise_p<0.05)==ncol(combinations)){
  # if all pairwise comparisons were significant
  nmds.plot<-nmds.plot+
    labs(caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else if(all.test.p>0.05 &sum(pairwise_p<0.05)==0){
  # if no comparisons were significant
  nmds.plot<-nmds.plot+
    labs(caption =  "No pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else if(all.test.p<0.05 &sum(pairwise_p<0.05)!=ncol(combinations)){
  # if some but not all comparisons were not significant
  nmds.plot<-nmds.plot+
    labs(caption =  paste("All pairwise comparisons except",
         gsub("_", " vs ", paste(names(which(pairwise_p>0.05)), collapse = ", ")),
         "were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons"))
}

nmds.plot.with.labels<-nmds.plot+
  ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE) # label each point by sample name


for(plot.type in seq_along(plot.types)){
  for(image.format in image.formats){
    ggsave(paste0("./images/diversity/nmds/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "ndms",dist.metric,
                        paste(custom.levels,collapse = '-'),
                        agglom.rank,truncationlvl,
                        sep = "-"),plot.types[plot.type],".",image.format),
           plot=get(paste0("nmds.",names(plot.types[plot.type]))),
           width = 4500,height = 3000,
           units = "px",dpi=300,device = image.format)
  }
}


## PCoA ####
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
  inner_join(.,custom.md, by="Sample")  

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
                     labels=unname(pretty.level.names),
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  labs(x=paste0("PCo 1 (", percent_exp[1],"%)"),
       y=paste0("PCo 2 (", percent_exp[2],"%)"),
       color="Host")+
  ggtitle(paste0("PCoA on ", beta.label, " between different rodents (",agglom.rank, " level)"))+
  mytheme # general theme

# Depending on the results of PERMANOVA, we should choose an appropriate
# label
if(all.test.p<0.05 &sum(pairwise_p<0.05)==ncol(combinations)){
  # if all pairwise comparisons were significant
  pcoa.plot<-pcoa.plot+
    labs(caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else if(all.test.p>0.05 &sum(pairwise_p<0.05)==0){
  # if no comparisons were significant
  pcoa.plot<-pcoa.plot+
    labs(caption =  "No pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else if(all.test.p<0.05 &sum(pairwise_p<0.05)!=ncol(combinations)){
  # if some but not all comparisons were not significant
  pcoa.plot<-pcoa.plot+
    labs(caption =  paste("All pairwise comparisons except",
                          gsub("_", " vs ", paste(names(which(pairwise_p>0.05)), collapse = ", ")),
                          "were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons"))
}

pcoa.plot.with.labels<-pcoa.plot+
  ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE) # label each point by sample name


for(plot.type in seq_along(plot.types)){
  for(image.format in image.formats){
    ggsave(paste0("./images/diversity/pcoa/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "pcoa",dist.metric,
                        paste(custom.levels,collapse = '-'),
                        agglom.rank,truncationlvl,
                        sep = "-"),plot.types[plot.type],".",image.format),
           plot=get(paste0("pcoa.",names(plot.types[plot.type]))),
           width = 4500,height = 3000,
           units = "px",dpi=300,device = image.format)
  }
}


save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  authorname,read.end.type,"beta-diversity",
  truncationlvl,agglom.rank,
  "workspace.RData",sep = "-")))

# PCA ####
# import rarefied dataframe
ps.q.df.pca.input<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      ps.q.df.pca.input.date_time,
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
           header = T)
# Convert into wide format
ps.q.df.wide.pca<-ps.q.df.pca.input%>%
  dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()


# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide.pca)<-ps.q.df.wide.pca$Sample
ps.q.df.wide.pca<-ps.q.df.wide.pca[,-1] # prune after this command 
ps.q.df.wide.pca<-as.matrix(ps.q.df.wide.pca)
# The first crucial step is to transform the compositional data:
# Use a log-ratio transformation, typically the centered log-ratio (clr)
# or isometric log-ratio (ilr) transformation.
# The clr transformation is often preferred for interpretability, 
# while ilr is used for statistical properties.
# We will use robust CLR (RCLR) with natural log
ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,method="rclr",logbase = exp)

#### >if you want to exclude specific samples
# pruned.samples<-c("PVO_15","PVO_19")
# ps.q.pruned<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%pruned.samples),]
# ps.q.pruned<-ps.q.pruned[,which(colSums(ps.q.pruned)!=0)]
# ps.q.df.wide.pca.tfm<-decostand(ps.q.pruned,method="rclr",logbase = exp)
####<

# calculate principal components
pca.q<-prcomp(ps.q.df.wide.pca.tfm)
str(pca.q)
dim(pca.q$x)

# reverse the signs
pca.q$rotation<- -1*pca.q$rotation

# display principal components (loadings)
head(pca.q$rotation)

# reverse the signs of the scores
pca.q$x<- -1*pca.q$x

# display the first six scores
head(pca.q$x)

## PCA Biplot ####
# biplot(pca.q,scale = 0)

#calculate total variance explained by each principal component
pca.q$sdev^2 / sum(pca.q$sdev^2)

#calculate total variance explained by each principal component
var_explained = pca.q$sdev^2 / sum(pca.q$sdev^2)

#create scree plot
qplot(seq_along(1:nrow(ps.q.df.wide.pca.tfm)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


## PCA Plot ####
PC1<-pca.q$x[,1]
PC2<-pca.q$x[,2]
perc.var<-round(100*summary(pca.q)$importance[2,1:2],2)

PC1.df<-as.data.frame(PC1)%>%rownames_to_column("Sample")
PC2.df<-as.data.frame(PC2)%>%rownames_to_column("Sample")

pca.plot<-PC1.df%>%
  left_join(PC2.df)%>%
  left_join(custom.md)%>%
  ggplot(aes(x=PC1,y=PC2,color=class,fill=class)) +
  geom_point(size=2)+ 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"),
       color="Host")+
  theme_bw()+
  ggtitle(paste0("PCA between different rodents (",agglom.rank, " level)"))+
  scale_color_manual(breaks = custom.levels,
                     labels=unname(pretty.level.names),
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  guides(fill="none")+
  mytheme # general theme

# Adjust the plot and filename if the data is pruned
if("pruned.samples"%in%ls()){
  pruned.or.not<-"-pruned"
  pca.plot<-pca.plot+
    labs(caption = paste("Removed samples:",paste(pruned.samples,collapse = ', ')))
}else{
  pruned.or.not<-""
}

pca.plot.with.labels<-pca.plot+
  ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE) # label each point by sample name

# Loop over the graphic devices
for(plot.type in seq_along(plot.types)){
  for(device.type in image.formats){
    ggsave(paste0("./images/diversity/pca/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "pca",
                        paste(custom.levels,collapse = '-'),
                        agglom.rank,truncationlvl,sep = "-"),
                  plot.types[plot.type],pruned.or.not,".",device.type),
           plot=get(paste0("pca.",names(plot.types[plot.type]))),
           width = 4500,height = 3000,
           units = "px",dpi=300,device = device.type)
  }
}



# find ASVs that contribute to PCs
aload<-abs(pca.q$rotation)
head(sweep(aload,2,colSums(aload),"/"))
colSums(sweep(aload, 2, colSums(aload), "/"))
colSums(aload)
aload[1,1]/colSums(aload)[1]

pc.df<-as.data.frame(sweep(aload,2,colSums(aload),"/")[,1:2])
lapply(pc.df,max)
max.ind<-lapply(pc.df,which.max)
pc.df[max.ind$PC1,]
pc.df[max.ind$PC2,]

pc.df%>%rownames_to_column(var=agglom.rank)%>%
  arrange(-PC2)%>%
  left_join(ps.q.agg[,agglom.rank])%>%
  distinct()%>%
  View

pc.df%>%arrange(-PC1)%>%head
pc.df%>%arrange(-PC2)%>%head
