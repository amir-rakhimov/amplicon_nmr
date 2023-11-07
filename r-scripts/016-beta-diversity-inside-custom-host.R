library(eeptools)
library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)

truncationlvl<-"234"
# agglom.rank<-"Genus"
agglom.rank<-"OTU"
read.end.type<-"single"

host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
comparison<-"sex"
comparison<-"strain"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

gg.title.taxon<-ifelse(agglom.rank=="OTU","(ASV level)",
                       paste0("(",agglom.rank," level)"))


if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.agg%>%
    filter(class=="NMR",Abundance!=0)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=age_calc(birthday,units = "years"))%>%
    mutate(age=round(age))
  
  min_boundary <- floor(min(ps.q.df.preprocessed$age)/5) * 5
  max_boundary <- ceiling(max(ps.q.df.preprocessed$age)/5) * 5
  
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(age_group=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                          include.lowest = TRUE))
}else if(host=="mice"){
  # select mice and add age groups
  custom.levels<-c("B6mouse",
                   "MSMmouse",
                   "FVBNmouse")
  ps.q.df.preprocessed<-ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
  ps.q.df.preprocessed$age_group<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                              ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
}


if (comparison=="age"){
    pretty.facet.labels<-names(table(ps.q.df.preprocessed$age_group))
    names(pretty.facet.labels)<-names(table(ps.q.df.preprocessed$age_group))
    custom.levels<-names(pretty.facet.labels)
    gg.labs.name<-"Age group"
    gg.title.groups<-"age groups"
    
}else if (comparison=="sex"){
    pretty.facet.labels<-
      c("F" = "Females",
        "M" = "Males")
    custom.levels<-names(pretty.facet.labels)
    pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]
    gg.labs.name<-"Host sex"
    gg.title.groups<-"groups"
}else if(comparison=="strain"){
  pretty.facet.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
  pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]
  gg.labs.name<-"Strain"
  gg.title.groups<-"strains"
}

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
dist.metric<-"bray"
if(dist.metric=="bray"){
  beta.label<-"Bray-Curtis dissimilarities"
}else { # in case of robust.aitchison, we need to substitute the dot
  beta.label<-paste(str_to_title(gsub("\\.", " ", dist.metric)),"distances")
}

permut.num<-1000 # number of permutations for PERMANOVA

ps.sampledata<-ps.q.df.preprocessed%>%
  select(c("Sample","class", "sex","birthday","age_group"))%>%
  distinct() # metadata

ps.q.df <-ps.q.df.preprocessed%>%
  select(Sample,OTU,Abundance,class,Taxon)


# Beta diversity ####
# find the smallest sample size
min.n_seqs.all<-ps.q.df%>%
  select(Sample,OTU,Abundance)%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# convert the data frame into wide format
if (agglom.rank=="OTU"){
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
  # dist.tfm<-decostand(x = ps.q.df.wide,"clr",pseudocount=1)
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

meta.dist<-dist.df%>%inner_join(ps.sampledata,.,by="Sample") # distances with metadata

## PERMANOVA ####
if (comparison=="age"){
  all.test<-adonis2(dist~age_group, # distances explained by class
                    data=meta.dist, 
                    permutations = permut.num) # permutation number
  
}else if (comparison=="sex"){
  all.test<-adonis2(dist~sex, # distances explained by class
                    data=meta.dist, 
                    permutations = permut.num) # permutation number
}else if(comparison=="strain"){
  all.test<-adonis2(dist~class, # distances explained by class
                    data=meta.dist, 
                    permutations = permut.num) # permutation number
}

all.test$F[1]
all.test.p<-all.test$`Pr(>F)`[1]
all.test.p

## Making pairwise comparisons ####
# we found a significant p value, but which groups are different?
pairwise_p<-numeric()
pairwise_F<-numeric()

# pairwise tests
combinations<-combn(custom.levels,2)
for (i in 1:ncol(combinations)) {
  if (comparison=="age"){
    test.data<-meta.dist%>%
      filter(age_group==combinations[1,i] | age_group==combinations[2,i]) # extract data for two groups
    
  }else if (comparison=="sex"){
    test.data<-meta.dist%>%
      filter(sex==combinations[1,i] | sex==combinations[2,i]) # extract data for two groups
  }else if(comparison=="strain"){
    test.data<-meta.dist%>%
      filter(class==combinations[1,i] | class==combinations[2,i]) # extract data for two groups
  }
  
  test.distances<-test.data%>%# extract only distances
    select(all_of(.$Sample))%>% # extract distances (corresponding to Sample vector)
    as.dist()
  
  if (comparison=="age"){
    adonis.test<-adonis2(test.distances~age_group,
                         data=test.data,
                         permutations=permut.num)
    
  }else if (comparison=="sex"){
    adonis.test<-adonis2(test.distances~sex,
                         data=test.data,
                         permutations=permut.num)
  }else if(comparison=="strain"){
    adonis.test<-adonis2(test.distances~class,
                         data=test.data,
                         permutations=permut.num)
  }
  
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

if(comparison=="age"){
  nmds.plot<-ggplot(nmds.scores,
                    aes(x=NMDS1,y=NMDS2,color=age_group, fill=age_group))
}else if (comparison=="sex"){
  nmds.plot<-ggplot(nmds.scores,
                    aes(x=NMDS1,y=NMDS2,color=sex, fill=sex))
}else if(comparison=="strain"){
  nmds.plot<-ggplot(nmds.scores,
                    aes(x=NMDS1,y=NMDS2,color=class, fill=class))
}

nmds.plot<-nmds.plot+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  theme_bw()+
  guides(fill="none")+
  # ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE)+ # add labels to samples+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  labs(color=gg.labs.name)+
  ggtitle(paste("nMDS on", beta.label, "between different",as.character(host.class[host]),gg.title.groups,"\n",
                gg.title.taxon))+
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

# if (comparison=="age"){
#   nmds.plot<-nmds.plot+
#     labs(color="Age groups")
#   
#   if(agglom.rank=="OTU"){
#     nmds.plot<-nmds.plot+
#       ggtitle(paste0("nMDS on ", beta.label, " between different ",as.character(host.class[host])," age groups 
#                    (ASV level)"))
#   }else{
#     nmds.plot<-nmds.plot+
#       ggtitle(paste0("nMDS on ", beta.label, " between different ",as.character(host.class[host])," age groups 
#                    (",agglom.rank, " level)"))
#   }
# }else if (comparison=="sex"){
#   nmds.plot<-nmds.plot+
#     labs(color="Host sex")
#   
#   if(agglom.rank=="OTU"){
#     nmds.plot<-nmds.plot+
#       ggtitle(paste0("nMDS on ", beta.label, " between different ",as.character(host.class[host])," groups 
#                   (ASV level)"))
#   }else{
#     nmds.plot<-nmds.plot+
#       ggtitle(paste0("nMDS on ", beta.label, " between different ",as.character(host.class[host])," groups 
#                   (",agglom.rank, " level)"))
#   }
# }else if(comparison=="strain"){
#   nmds.plot<-nmds.plot+
#     labs(color="Strain")
#   if(agglom.rank=="OTU"){
#     nmds.plot<-nmds.plot+
#     ggtitle(paste0("nMDS on ", beta.label, " between different mouse strains 
#                    (ASV level)"))
#   }else{
#     nmds.plot<-nmds.plot+
#       ggtitle(paste0("nMDS on ", beta.label, " between different mouse strains 
#                    (",agglom.rank, " level)"))
#   }
# }

ggsave(paste0("./images/diversity/nmds/",
              paste(Sys.Date(),"ndms",dist.metric,host,
                    comparison,agglom.rank,truncationlvl,sep = "-"),".png"),
       plot=nmds.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/nmds/",
              paste(Sys.Date(),"ndms",dist.metric,host,
                    comparison,agglom.rank,truncationlvl,sep = "-"),".tiff"),
       plot=nmds.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")



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

if(comparison=="age"){
  pcoa.plot<-ggplot(pcoa.positions,
                    aes(x=pcoa1,y=pcoa2,color=age_group,
                        fill=age_group))
}else if (comparison=="sex"){
  pcoa.plot<-ggplot(pcoa.positions,
                    aes(x=pcoa1,y=pcoa2,color=sex,
                        fill=sex))
}else if(comparison=="strain"){
  pcoa.plot<-ggplot(pcoa.positions,
                    aes(x=pcoa1,y=pcoa2,color=class,fill=class))
}

pcoa.plot<-pcoa.plot+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  theme_bw()+
  guides(fill="none")+
  labs(x=paste0("PCo 1 (", percent_exp[1],"%)"),
       y=paste0("PCo 2 (", percent_exp[2],"%)"),
       color=gg.labs.name
       # ,
       # caption = "All pairwise comparisons were significant using<br>ADONIS at p=0.05 using Benjamini-Hochberg correction for multiple comparisons"
  )+
  ggtitle(paste("PCoA on", beta.label, "between different",as.character(host.class[host]),gg.title.groups, "\n",
                gg.title.taxon))+
  # ggrepel::geom_text_repel(aes(label=Sample),show.legend = FALSE)+ # add labels to samples+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
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

# if (comparison=="age"){
#   pcoa.plot<-pcoa.plot+
#     labs(color="Age groups")
# 
#   if(agglom.rank=="OTU"){
#     pcoa.plot<-pcoa.plot+
#       ggtitle(paste0("PCoA on ", beta.label, " between different naked mole-rat age groups
#                    (ASV level)"))
#   }else{
#     pcoa.plot<-pcoa.plot+
#       ggtitle(paste0("PCoA on ", beta.label, " between different naked mole-rat age groups
#                    (",agglom.rank, " level)"))
#   }
# }else if (comparison=="sex"){
#   pcoa.plot<-pcoa.plot+
#     labs(x=paste0("PCo 1 (", percent_exp[1],"%)"),
#          y=paste0("PCo 2 (", percent_exp[2],"%)"),
#          color="Host sex")
# 
#   if(agglom.rank=="OTU"){
#     pcoa.plot<-pcoa.plot+
#       ggtitle(paste0("PCoA on ", beta.label, " between different naked mole-rat groups
#                    (ASV level)"))
#   }else{
#     pcoa.plot<-pcoa.plot+
#       ggtitle(paste0("PCoA on ", beta.label, " between different naked mole-rat groups
#                    (",agglom.rank, " level)"))
#   }
# }else if(comparison=="strain"){
#   if(agglom.rank=="OTU"){
#     pcoa.plot<-pcoa.plot+
#       ggtitle(paste0("PCoA on ", beta.label, " between different mouse strains
#                    (ASV level)"))+
#       labs(color="Strain")
#   }else{
#     pcoa.plot<-pcoa.plot+
#       ggtitle(paste0("PCoA on ", beta.label, " between different mouse strains
#                    (",agglom.rank, " level)"))+
#       labs(color="Strain")
#   }
# }

ggsave(paste0("./images/diversity/pcoa/",
              paste(Sys.Date(),"pcoa",dist.metric,host,
                    comparison,agglom.rank,truncationlvl,sep = "-"),".png"),
       plot=pcoa.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")

ggsave(paste0("./images/diversity/pcoa/",
              paste(Sys.Date(),"pcoa",dist.metric,host,
                    comparison,agglom.rank,truncationlvl,sep = "-"),".tiff"),
       plot=pcoa.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")

# PCA ####
ps.q.df.wide.centered<-scale(ps.q.df.wide,scale=F,center=T)
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
biplot(pca.q,scale = 0)

#calculate total variance explained by each principal component
pca.q$sdev^2 / sum(pca.q$sdev^2)

#calculate total variance explained by each principal component
var_explained = pca.q$sdev^2 / sum(pca.q$sdev^2)

#create scree plot
qplot(seq_along(1:nrow(ps.q.df.wide)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)



PC1<-pca.q$x[,1]
PC2<-pca.q$x[,2]
perc.var<-round(100*summary(pca.q)$importance[2,1:2],2)

if(comparison=="age"){
  pca.plot<-ggplot(ps.sampledata,
                    aes(x=PC1,y=PC2,color=age_group,
                        fill=age_group))
}else if (comparison=="sex"){
  pca.plot<-ggplot(ps.sampledata,
                    aes(x=PC1,y=PC2,color=sex,
                        fill=sex))
}else if(comparison=="strain"){
  pca.plot<-ggplot(ps.sampledata,
                    aes(x=PC1,y=PC2,color=class,fill=class))
}

pca.plot<-pca.plot +
  geom_point(size=2)+ 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"),
       color=gg.labs.name)+
  theme_bw()+
  ggtitle(paste("PCA between different",as.character(host.class[host]),gg.title.groups, "\n",
                gg.title.taxon))+
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
ggsave(paste0("./images/diversity/pca/",
              paste(Sys.Date(),"pca",host,
                    comparison,agglom.rank,truncationlvl,sep = "-"),".png"),
       plot=pca.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "png")

ggsave(paste0("./images/diversity/pca/",
              paste(Sys.Date(),"pca",host,
                    comparison,agglom.rank,truncationlvl,sep = "-"),".tiff"),
       plot=pca.plot,
       width = 4500,height = 3000,
       units = "px",dpi=300,device = "tiff")


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
