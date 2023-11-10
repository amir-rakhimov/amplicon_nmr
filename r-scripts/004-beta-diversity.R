library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"

read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))

pretty.facet.labels<-
  c("NMR" = "*Heterocephalus glaber*", # better labels for facets
    "B6mouse" = "B6 mouse",
    # "MSMmouse" = "MSM/Ms mouse",
    # "FVBNmouse" = "FVB/N mouse",
    "DMR" = "*Fukomys Damarensis*",
    "hare" = "*Lepus europaeus*",
    "rabbit" = "*Oryctolagus cuniculus*",
    "spalax" = "*Nannospalax leucodon*",
    "pvo" = "*Pteromys volans orii*"
    # "NMRwt"="Wild *Heterocephalus glaber*",
)

custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
# Using Polychrome package
set.seed(1)
# custom.colors<- createPalette(length(custom.levels),
#               seedcolors = rainbow(6))
custom.colors<- 
  createPalette(length(custom.levels),
                seedcolors = c("#B22222", "#0000FF","#006400", 
                               "#FF8C00","#5D478B", "#00FFFF"))
custom.colors<-unname(custom.colors)
swatch(custom.colors)
# custom colors for scale_fill_manual (maybe not needed)
# custom.fill<-c("dodgerblue", "seagreen","indianred")

# custom labels for scale_color_manual
custom.color.labels<-unname(pretty.facet.labels)

# Choose distance metric
dist.metric<-"robust.aitchison"
dist.metric<-"jaccard"
dist.metric<-"canberra"
dist.metric<-"uni.unw"
dist.metric<-"bray"
if(dist.metric=="bray"){
  beta.label<-"Bray-Curtis dissimilarities"
}else if(dist.metric=="uni.unw"){
  beta.label<-"Unweighted UniFrac distances"
}else { # in case of robust.aitchison, we need to substitute the dot
  beta.label<-paste(str_to_title(gsub("\\.", " ", dist.metric)),"distances")
}

permut.num<-1000 # number of permutations for PERMANOVA

ps.sampledata<-ps.q.agg%>%
  select(c("Sample","class","sex","birthday"))%>%
  distinct() # metadata

ps.q.df <-ps.q.agg%>%
  # filter(Taxon.bp!="Remainder (Mean abundance < 1%)")%>%
  select(Sample,OTU,Abundance,class,Taxon)

# ps.q.df <-ps.q.agg.rel%>%
#   select(Sample,OTU,Abundance,class,Taxon)
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

# convert the data frame into wide format
if(agglom.rank=="OTU"){
  # convert the data frame into wide format
  ps.q.df.wide<-ps.q.df%>%
    select(-Taxon)%>%
    pivot_wider(names_from = "OTU", # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()
}else{
  ps.q.df.wide<-ps.q.df%>%
    select(-OTU)%>%
    pivot_wider(names_from = "Taxon", # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()
  
}


# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  
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
                # ,
                # pseudocount=1
                )
}else if(dist.metric=="uni.unw"){ # TODO: fix unifrac
  dist<-UniFrac(ps.q, weighted = FALSE)
}

dist.df<-dist%>% # convert distances into tibble
  as.matrix()%>%
  as_tibble()
dist.df<-cbind(labels(dist),dist.df) # create a column of sample ids
colnames(dist.df)[1]<-"Sample"

meta.dist<-dist.df%>%inner_join(ps.sampledata,.,by="Sample") # distances with metadata

## PERMANOVA ####
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
  labs(color="Host")+
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

# Depending on the results of PERMANOVA, we should choose an appropriate
# label
if(all.test.p<0.05){
  # if all pairwise comparisons were significant
  nmds.plot<-nmds.plot+
    labs(caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else if(all.test>0.05 &length(pairwise_p)==0){
  # if no comparisons were significant
  nmds.plot<-nmds.plot+
    labs(caption =  "No pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else{
  # if some but not all comparisons were not significant
  nmds.plot<-nmds.plot+
    labs(caption =  paste("All pairwise comparisons except",
         gsub("_", " vs ", paste(names(which(pairwise_p<0.001)), collapse = ", ")),
         "were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons"))
}

ggsave(paste0("./images/diversity/nmds/",
              paste(Sys.Date(),"ndms",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,sep = "-"),".png"),
       plot=nmds.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")

ggsave(paste0("./images/diversity/nmds/",
              paste(Sys.Date(),"ndms",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,sep = "-"),".tiff"),
       plot=nmds.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")

nmds.plot.with.labels<-nmds.plot+
  ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE) # label each point by sample name
ggsave(paste0("./images/diversity/nmds/",
              paste(Sys.Date(),"ndms",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,"with-labels",sep = "-"),".png"),
       plot=nmds.plot.with.labels,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/nmds/",
              paste(Sys.Date(),"ndms",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,"with-labels",sep = "-"),".tiff"),
       plot=nmds.plot.with.labels,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")

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
       color="Host")+
  ggtitle(paste0("PCoA on ", beta.label, " between different rodents (",agglom.rank, " level)"))+
  theme(plot.title = element_text(size = 27),
    axis.text.x = element_text(angle=0,size=20,hjust=1),
    axis.text.y = element_text(size=20),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.position = "right",
    legend.text = ggtext::element_markdown(size = 20),
    plot.caption = ggtext::element_markdown(hjust = 0, size=20),
    plot.caption.position = "plot")

# Depending on the results of PERMANOVA, we should choose an appropriate
# label
if(all.test.p<0.05){
  # if all pairwise comparisons were significant
  pcoa.plot<-pcoa.plot+
    labs(caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else if(all.test>0.05 &length(pairwise_p)==0){
  # if no comparisons were significant
  pcoa.plot<-pcoa.plot+
    labs(caption =  "No pairwise comparisons were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons")
}else{
  # if some but not all comparisons were not significant
  pcoa.plot<-pcoa.plot+
    labs(caption =  paste("All pairwise comparisons except",
                          gsub("_", " vs ", paste(names(which(pairwise_p<0.001)), collapse = ", ")),
                          "were significant using<br>ADONIS at p=0.05 with Benjamini-Hochberg correction for multiple comparisons"))
}

ggsave(paste0("./images/diversity/pcoa/",
              paste(Sys.Date(),"pcoa",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,sep = "-"),".png"),
       plot=pcoa.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/pcoa/",
              paste(Sys.Date(),"pcoa",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,sep = "-"),".tiff"),
       plot=pcoa.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")

pcoa.plot.with.labels<-pcoa.plot+
  ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE) # label each point by sample name

ggsave(paste0("./images/diversity/pcoa/",
              paste(Sys.Date(),"pcoa",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,"with-labels",sep = "-"),".png"),
       plot=pcoa.plot.with.labels,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/pcoa/",
              paste(Sys.Date(),"pcoa",dist.metric,
                    paste(custom.levels,collapse = '-'),
                    agglom.rank,truncationlvl,"with-labels",sep = "-"),".tiff"),
       plot=pcoa.plot.with.labels,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")


# PCA ####
# import rarefied dataframe
ps.q.df.pca.input<-read.table(paste0("./rtables/",authorname,"/ps.q.df.",
                  rare.status,".",filter.status,"-",agglom.rank,"-",
                  paste(custom.levels,collapse = '-'),".tsv"),
           header = T)
# Convert into wide format
if(agglom.rank=="OTU"){
  ps.q.df.wide<-ps.q.df.pca.input%>%
    select(-class,-sex,-birthday)%>%
    pivot_wider(names_from = "OTU", # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()
}else{
  ps.q.df.wide<-ps.q.df.pca.input%>%
    select(-class,-sex,-birthday)%>%
    pivot_wider(names_from = "Taxon", # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()
}

# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-1] # prune after this command 
ps.q.df.wide<-as.matrix(ps.q.df.wide)
ps.q.df.wide.centered<-scale(ps.q.df.wide,scale=F,center=T)

#### >if you want to exclude specific samples
pruned.samples<-c("PVO_15","PVO_19")
ps.q.pruned<-ps.q.df.wide[-which(rownames(ps.q.df.wide)%in%pruned.samples),]
ps.q.pruned<-ps.q.pruned[,which(colSums(ps.q.pruned)!=0)]
ps.q.df.wide.centered<-scale(ps.q.pruned,scale=F,center=T)
####<

ps.q.df.wide.centered.scaled<-scale(ps.q.df.wide.centered,scale=T,center=F)
# calculate principal components
pca.q<-prcomp(ps.q.df.wide.centered.scaled)
str(pca.q)
dim(pca.q$x)

# reverse the signs
pca.q$rotation<- -1*pca.q$rotation

# display principal components (loadings)
head(pca.q$rotation)

# reverse th signs of the scores
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
qplot(seq_along(1:nrow(ps.q.df.wide.centered.scaled)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)


## PCA Plot ####
PC1<-pca.q$x[,1]
PC2<-pca.q$x[,2]
perc.var<-round(100*summary(pca.q)$importance[2,1:2],2)

pca.plot<-ggplot(ps.sampledata[ps.sampledata$Sample%in%
                                 rownames(ps.q.df.wide.centered.scaled),],
                 aes(x=PC1,y=PC2,color=class,fill=class)) +
  geom_point(size=2)+ 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"),
       color="Host")+
  theme_bw()+
  ggtitle(paste("PCA between different rodents (",agglom.rank, " level)"))+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  guides(fill="none")+
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

if(length(pruned.samples)!=0){
  pca.plot<-pca.plot+
    labs(caption = paste("Removed samples:",paste(pruned.samples,collapse = ', ')))
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,"pruned",sep = "-"),".png"),
         plot=pca.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "png")
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,"pruned",sep = "-"),".tiff"),
         plot=pca.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "tiff")
}else{
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,sep = "-"),".png"),
         plot=pca.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "png")
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,sep = "-"),".tiff"),
         plot=pca.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "tiff")
}


pca.plot.with.labels<-pca.plot+
  ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE) # label each point by sample name

if(length(pruned.samples)!=0){
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,"with-labels-pruned",sep = "-"),".png"),
         plot=pca.plot.with.labels,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "png")
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,"with-labels-pruned",sep = "-"),".tiff"),
         plot=pca.plot.with.labels,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "tiff")
}else{
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,"with-labels",sep = "-"),".png"),
         plot=pca.plot.with.labels,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "png")
  ggsave(paste0("./images/diversity/pca/",
                paste(Sys.Date(),"pca",
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,"with-labels",sep = "-"),".tiff"),
         plot=pca.plot.with.labels,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = "tiff")
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

pc.df%>%rownames_to_column(var="OTU")%>%
  arrange(-PC2)%>%
  left_join(ps.q.agg[,c("OTU","Taxon")])%>%
  distinct()%>%
  View


