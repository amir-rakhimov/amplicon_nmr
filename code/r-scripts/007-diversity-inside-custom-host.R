# Alpha and beta diversity inside custom host #### 
library(vegan)
library(tidyverse)
library(Polychrome)
# Import data ####
ps.q.df.preprocessed.date_time<-"20240524_13_58_11"
# 20240426_22_00_04 preprocessed df for all hosts, genus level
# "20240524_13_58_11" preprocessed df for NMR, OTU level
# 20240809_13_18_49 preprocessed df for NMR, Genus level
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"OTU"
read.end.type<-"single"
rare.status<-"rare"
rdafiles.directory<-"./output/rdafiles/use"
rtables.directory<-file.path("./output/rtables",authorname)

ps.q.agg.date_time<-"20240620_12_38_18"
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))

# 20240620_12_38_18 ps.q.agg.date_time OTU level, all hosts
# 20240620_12_40_41 ps.q.agg.date_time genus level, all hosts
# 20240917_10_54_12 ps.q.agg.date_time family level, all hosts
# 20240917_21_29_36 ps.q.agg.date_time phylum level, all hosts
# Import metadata
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
image.formats<-c("png","tiff")
filter.status<-"nonfiltered"
# If you compare NMRs or mice (inside host)
host<-"NMR"
# host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"

gg.title.taxon<-ifelse(agglom.rank=="OTU","(ASV level)",
                       paste0("(",agglom.rank," level)"))

if(host=="NMR"){
  host.labels<-
    c("NMR" = "*Heterocephalus glaber*")
}else if(host=="mice"){
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}


# load the output of 003-phyloseq-rarefaction-filtering.R file
ps.q.df.preprocessed<-read.table(
  file.path(rtables.directory,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
      paste(names(host.labels),collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

if(host=="NMR"){
  # Load a table with ages
  custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))
}else if(host=="mice"){
  # this is an if/else statement to create the agegroup column
  custom.md%>%
    filter(class%in%c("B6mouse","FVBNmouse","MSMmouse"))%>%
    mutate(agegroup=case_when(
      class == "B6mouse" ~ "B6",
      class == "MSMmouse" & grepl("2020", birthday) ~ "MSM/Ms old",
      class == "MSMmouse" & as.Date(birthday) > as.Date("2021-01-01") ~ "MSM/Ms young",
      class == "FVBNmouse" ~ "FVB/N young",
      TRUE ~ NA # NA if none of the conditions are satisfied
    ))
}


if (comparison=="age"){
  pretty.level.names<-names(table(custom.md$old_agegroup))
  names(pretty.level.names)<-names(table(custom.md$agegroup))
  custom.levels<-names(pretty.level.names)
  gg.labs.name<-"Age group"
  gg.title.groups<-"age groups"
  
}else if (comparison=="sex"){
  pretty.level.names<-
    c("F" = "Females",
      "M" = "Males")
  custom.levels<-names(pretty.level.names)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  gg.labs.name<-"Host sex"
  gg.title.groups<-"groups"
}else if(comparison=="strain"){
  pretty.level.names<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.level.names),custom.md$class)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  gg.labs.name<-"Strain"
  gg.title.groups<-"strains"
}else if(comparison =="colony"){
  relations<-read.table("./data/metadata/pooled-metadata/nmr-relations.tsv",
                                   header = T,
                                   sep = "\t")
  custom.md<-custom.md%>%
    left_join(relations,by="Sample")
  pretty.level.names<-names(table(custom.md$colony))
  names(pretty.level.names)<-names(table(custom.md$colony))
  custom.levels<-names(pretty.level.names)
  gg.labs.name<-"Colony"
  gg.title.groups<-"colonies"
}


metric.labs=c('sobs'="Richness (Observed species)",
              'shannon' = "Shannon",
              # 'simpson' = "Simpson",
              'invsimpson' = "Inverse Simpson")
plot.metrics<-c("sobs","shannon", # "simpson",
                "invsimpson") # metrics to plot
div.indices<-c("sobs","shannon",# "simpson",
               "invsimpson")

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

# Theme for alpha diversity plots
alpha.plot.theme<-theme(plot.margin=unit(c(1,1,1,2), 'cm'),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.text.x = ggtext::element_markdown(angle=45,hjust=1,size=18),
                                      axis.text.y = element_text(size=20),
                                      axis.title = element_text(size = 20),
                                      strip.text.x = element_text(size=20),
                                      plot.title = element_text(size = 27),
                                      legend.text = element_text(size = 20),
                                      legend.title = element_text(size = 25))
# Theme for scatter plots (PCA, nMDS, PCoA)
scatter.plot.theme<-theme(
  axis.text.x = element_text(angle=0,size=20,hjust=1),
  axis.text.y = element_text(size=20),
  axis.title = element_text(size = 20),
  plot.title = element_text(size = 27),
  legend.text = ggtext::element_markdown(size = 20),
  legend.title = element_text(size = 25),
  legend.position = "right",
  plot.caption = ggtext::element_markdown(hjust = 0, size=20),
  plot.caption.position = "plot")

# Alpha diversity ####
## Compute alpha diversity metrics ####
all.div<-ps.q.df.preprocessed%>%
  group_by(Sample)%>%
  reframe(sobs=specnumber(Abundance), # richness (num of species)
          shannon=diversity(Abundance,index = "shannon"),
          # simpson=diversity(Abundance, index="simpson"),
          invsimpson=diversity(Abundance, index="invsimpson"), # inverse simpson
          tot=sum(Abundance))%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()


## Alpha diversity tests ####
kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
colnames(kt.results)<-div.indices
rownames(kt.results)<-c("statistic", "pvalue")

combinations<-combn(custom.levels,2) # all unique pairwise combinations
w.results<-data.frame(matrix(nrow = ncol(combinations),ncol=length(div.indices))) # ncol(combinations) pairwise comparisons
colnames(w.results)<-div.indices
w.results<-array(dim = c(length(table(combinations[1,])), # nrows
                         length(table(combinations[2,])), # ncols
                         length(div.indices)), # num of 2D arrays (stacking) 
                 dimnames = list(NULL, NULL, div.indices))



for (div.metric in div.indices) {
  metric.ds<-all.div%>%
    left_join(custom.md)%>%
    filter(metric==div.metric)%>%
    distinct()
  # perform kruskal-wallis test
  if(comparison=="age"){
    kt<-kruskal.test(value~agegroup,data=metric.ds)
    
  }else if (comparison=="sex"){
    kt<-kruskal.test(value~sex,data=metric.ds)
    
  }else if(comparison=="strain"){
    kt<-kruskal.test(value~class,data=metric.ds)
    
  }
  
  kt.results["statistic",div.metric]<- kt$statistic
  kt.results["pvalue",div.metric]<- kt$p.value
  # low pvalue indicates statistical difference between some groups
  # But which groups are different from others?
  # Perform pairwise wilcoxon test
  
  # NB: You must not run pairwise wilcoxon test unless kruskal wallis results 
  # are significant
  if(kt$p.value<0.05){
    if(comparison=="age"){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds$agegroup,
                                   p.adjust.method = "BH",
                                   exact=FALSE)
    }else if (comparison=="sex"){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds$sex,
                                   p.adjust.method = "BH",
                                   exact=FALSE)
      
    }else if(comparison=="strain"){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds$class,
                                   p.adjust.method = "BH",
                                   exact=FALSE)
    }
    
    w.results[,,div.metric]<-w.test$p.value
    
  }else(
    w.results[,,div.metric]<-matrix(data = "n.s.",
                                    nrow = nrow(w.test$p.value),
                                    ncol = ncol(w.test$p.value))
    
  )
  dimnames(w.results)[[1]]<-dimnames(w.test$p.value)[[1]] # change rownames of w.results
  dimnames(w.results)[[2]]<-dimnames(w.test$p.value)[[2]] # change colnames of w.results
}

kt.results
w.results
stopifnot(all(kt.results[2,]<0.05))

# convert multidimensional array into data frame
# w.results.df<-w.results%>%
#   as.data.frame.table()%>%
#   drop_na()
# colnames(w.results.df)<-c("class1","class2","div.metric","p.value")


max.values<-all.div%>%
  group_by(metric)%>%
  summarise(max_val=max(value))

# Prepare data for plotting
if(comparison=="age"){
  custom.md$agegroup<-factor(custom.md$agegroup,levels=custom.levels)
}else if (comparison=="sex"){
  custom.md$sex<-factor(custom.md$sex,levels=custom.levels)
  
}else if(comparison=="strain"){
  custom.md$agegroup<-factor(custom.md$class,levels=custom.levels)
}else if(comparison=="colony"){
  custom.md$colony<-factor(custom.md$colony,levels=custom.levels)
}

## Plot alpha diversity metrics ####
if(comparison=="age"){
  div.plot<-all.div%>%
    filter(metric%in%plot.metrics)%>%
    left_join(custom.md)%>%
    ggplot(aes(x=factor(agegroup, level=custom.levels),
               y=value,
               fill=factor(agegroup)))
}else if (comparison=="sex"){
  div.plot<-all.div%>%
    filter(metric%in%plot.metrics)%>%
    left_join(custom.md)%>%
    ggplot(aes(x=factor(sex, level=custom.levels),
               y=value,
               fill=factor(sex)))
}else if(comparison=="strain"){
  div.plot<-all.div%>%
    filter(metric%in%plot.metrics)%>%
    left_join(custom.md)%>%
    ggplot(aes(x=factor(class,level=custom.levels),
               y=value,
               fill=factor(class)))
}else if(comparison == "colony"){
  div.plot<-all.div%>%
    filter(metric%in%plot.metrics)%>%
    left_join(custom.md)%>%
    ggplot(aes(x=reorder(colony,-value),y=value,fill=factor(colony)))
                   # aes(x=factor(class,
                   #              level=custom.levels),y=value,fill=factor(colony)))
}



div.plot<-div.plot+
  geom_boxplot(show.legend = FALSE,outliers = FALSE)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             ncol=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+ 
  labs(color=gg.labs.name)+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=custom.levels)+ # rename boxplot labels (x axis)
  scale_fill_manual(values = custom.fill)+
  alpha.plot.theme + # theme for alpha diversity plots
  theme(legend.position = "none")+
  ggtitle(paste("Alpha diversity of the gut microbiota of different",as.character(host.class[host]),gg.title.groups,
                gg.title.taxon))

div.plot.jitter<-div.plot+
  geom_jitter(width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=3,
              show.legend = FALSE)

div.plot.dots<-div.plot+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) # add dots

for(image.format in image.formats){
  ggsave(paste0("./images/diversity/alpha/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "alpha-jitter",
                      paste(plot.metrics,collapse = "-"),
                      host,comparison,agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.jitter,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = image.format)
  ggsave(paste0("./images/diversity/alpha/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "alpha-dots",
                      paste(plot.metrics,collapse = "-"),
                      host,comparison,agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.dots,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = image.format)
}

# Per-sample plot for ages ####
## Setup sample levels ####
if(comparison=="age"){
  sample.levels<-custom.md%>%
    ungroup()%>%
    select(Sample,age)%>%
    arrange(age)%>%
    distinct()%>%
    mutate(NewSample=paste0(Sample," (",age,")"))%>%
    mutate(Sample=factor(Sample,levels=Sample),
           NewSample=factor(NewSample,levels=NewSample))
  # You can plot diversity metrics per sample or average metrics per age.
  # When plotting average metrics per age, change the x axis text angle from 45 to 0.
  per.sample.div.plot<-all.div%>%
    left_join(sample.levels)%>%
    group_by(age,metric)%>%
    mutate(mean_metric=mean(value))%>% 
    ggplot(aes(x=age,y=mean_metric,colour=age))+ # when plotting average metric per age
    # ggplot(aes(x=NewSample,y=value,colour=age))+
    geom_point(size=4,stat = "identity")+
    # facet_wrap(~factor(metric, # reorder facets
    #                    levels=plot.metrics),
    #            nrow=length(plot.metrics),
    #            scales="free_y", # free y axis range
    #            labeller = as_labeller(metric.labs))+ # rename facets
    # Below is facets when we plot average values per age
    facet_wrap(~factor(metric, # reorder facets
                       levels=plot.metrics),
               nrow=length(plot.metrics),
               scales="free_y", # free y axis range
               labeller = as_labeller(sapply(metric.labs, function(x) paste("Average",x))))+ # rename facets
    scale_color_viridis_c(option="D")+
    scale_x_continuous(breaks = seq(min(sample.levels$age), max(sample.levels$age), by = 1))+
    theme_bw()+
    alpha.plot.theme + # theme for alpha diversity plots
    theme(legend.position = "right")+
    ggtitle(paste("Alpha diversity of the naked mole-rat gut microbiota across age",gg.title.taxon))
  
  for(image.format in image.formats){
    ggsave(paste0("./images/diversity/alpha/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "alpha-per-sample",
                        paste(plot.metrics,collapse = "-"),
                        host,comparison,agglom.rank,truncationlvl,
                        sep = "-"),".",image.format),
           plot=per.sample.div.plot,
           width = 6000,height = 3000,
           units = "px",dpi=300,device = image.format)
  }
}
  
# Add significance stars if we have significant results ####
if(length(which(as.vector(w.results)<0.05))>0){
  coord.combs<-combn(seq_along(custom.levels),2)
  
  # First, find significant results (stars)
  stars.list<-matrix(NA,nrow=length(div.indices),ncol = ncol(coord.combs))
  for (k in 1:dim(w.results)[3]) { # loop over each sub array
    for (i in 1:dim(w.results)[1]) { # then each row
      for (j in 1:dim(w.results)[2]) { # then each column
        if (!is.na(w.results[i,j,k])) { # extract non-NA values
          ith.row<-rownames(w.results[,,k])[i] 
          jth.col<-colnames(w.results[,,k])[j]
          w.val<- w.results[i,j,k]
          
          # e.g. ith.row="NMR"
          # jth.col="hare"
          # find in custom.levels the positions that correspond to c("NMR","hare")
          # the result: level.position = c(1,4)
          levels.position<-which(custom.levels%in%c(ith.row,jth.col))
          # find these positions in the coord.combs to get the column 
          # that stores levels.position
          # result: level.col=3 
          level.col<-which(apply(coord.combs,2,function(x) 
            return(all(x==levels.position))))
          # assign the w.val to stars.list:kth row is diversity metric
          # level.col column is the combination of levels
          # we must assign either "*" or "n.s." depeding on w.result[i,j,k]
          stars.list[k,level.col]<-ifelse(w.val!="n.s.",ifelse(w.val<0.05,"*","n.s."),w.val)
        }
      }
    }
  }
  
  stars.list<-as.vector(t(stars.list))
  
  stars<-tibble(
    metric=factor(rep(plot.metrics,each=ncol(coord.combs)), # each plot metric repeated by the number of combinations
                  levels = plot.metrics), # names of our metrics
    label=stars.list
  )
  star.indices<-which(stars$label=="*") # only significant results
  
  freqs<-as.data.frame(table(stars))%>%filter(label=="*")
  yvalues<-c()
  for (i in seq_along(freqs)){
    vec<-seq(from=1, by=0.08,length.out=freqs[i,"Freq"])*
      as.numeric(max.values[max.values$metric==freqs[i,"metric"],"max_val"])
    yvalues<-c(yvalues,vec)
  }
  
  
  # the horizontal lines
  # merge two vectors: two rows
  # c() turns them into a vector
  
  # y coordinates
  # multiply vectors to get a matrix: nrows is the number of comparisons
  # ncols is the number of plot.metrics
  # seq is the sequence of values that will multiply maxvalues
  
  # start and end values can be found from pairwise combinations
  # x values have dimensions: num of metrics * num of pairwise comparisons
  xvalues<-rep(coord.combs[1,],
               length=length(div.indices)*ncol(coord.combs)) # first row is x start values
  xendvalues<-rep(coord.combs[2,],
                  length=length(div.indices)*ncol(coord.combs)) # second row is x end values
  
  # select only significant results
  xvalues<-xvalues[star.indices]
  xendvalues<-xendvalues[star.indices]
  
  horizontal.lines<-tibble(
    metric=factor(rep(plot.metrics,each=ncol(coord.combs)),
                  levels = plot.metrics)[star.indices], # names of our metrics
    x = xvalues,       #1 2 1 1 2 1 1 2 1 1 2 1
    xend = xendvalues, #2 3 3 2 3 3 2 3 3 2 3 3
    y =yvalues,
    yend = yvalues
  )
  
  
  
  
  # labels depend on our tests (statistical significance)
  xstars<-(xvalues+xendvalues)/2
  stars<-stars%>% filter(label=="*")%>%
    mutate(x= xstars,
           y =c(yvalues)*1.01)
  
  
  
  newplot<-div.plot+
    geom_segment(data=horizontal.lines, # add horizontal.lines of significance
                 aes(x=x, xend=xend, y=y, yend=yend),
                 inherit.aes = FALSE)+# no conflict with different fills
    geom_text(data = stars,aes(x=x, y=y, label=label),
              inherit.aes = FALSE,size=10) # add stars 
  
  for(image.format in image.formats){
    ggsave(paste0("./images/diversity/alpha/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "alpha",
                        paste(plot.metrics,collapse = "-"),
                        host,comparison,agglom.rank,truncationlvl,"signif",
                        sep = "-"),".",image.format),
           plot=newplot,
           width = 6000,height = 3000,
           units = "px",dpi=300,device = image.format)
  }
}



# Beta diversity ####
set.seed(1)
custom.colors<- 
  createPalette(length(custom.levels),
                seedcolors = c("#B22222", "#0000FF","#006400", 
                               "#FF8C00","#5D478B", "#00FFFF"))
custom.colors<-unname(custom.colors)
swatch(custom.colors)

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

# find the smallest sample size
min.n_seqs.all<-ps.q.df.preprocessed%>%
  select(Sample,Abundance,class,all_of(agglom.rank))%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# convert the data frame into wide format
ps.q.df.wide<-ps.q.df.preprocessed%>%
  select(Sample,Abundance,class,all_of(agglom.rank))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()


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

meta.dist<-dist.df%>%
  inner_join(custom.md,.,by="Sample")%>% # distances with metadata
  ungroup()

## PERMANOVA ####
if (comparison=="age"){
  all.test<-adonis2(dist~agegroup, # distances explained by class
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
      filter(agegroup==combinations[1,i] | agegroup==combinations[2,i]) # extract data for two groups
    
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
    adonis.test<-adonis2(test.distances~agegroup,
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
  inner_join(.,custom.md, by="Sample") 

if(comparison=="age"){
  nmds.plot<-ggplot(nmds.scores,
                    aes(x=NMDS1,y=NMDS2,color=agegroup, fill=agegroup))
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
                     labels=unname(pretty.level.names),
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  labs(color=gg.labs.name)+
  ggtitle(paste("nMDS on", beta.label, "between different",as.character(host.class[host]),gg.title.groups,"\n",
                gg.title.taxon))+
  scatter.plot.theme  # theme for scatter plots (PCA, PCoA, nMDS)



for(image.format in image.formats){
  ggsave(paste0("./images/diversity/nmds/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "ndms",
                      dist.metric,host,
                      comparison,agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=nmds.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = image.format)
}


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
  inner_join(.,custom.md, by="Sample")  

if(comparison=="age"){
  pcoa.plot<-ggplot(pcoa.positions,
                    aes(x=pcoa1,y=pcoa2,color=agegroup,
                        fill=agegroup))
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
                     labels=unname(pretty.level.names),
                     values = custom.colors)+
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  scatter.plot.theme  # theme for scatter plots (PCA, PCoA, nMDS)

for(image.format in image.formats){
  ggsave(paste0("./images/diversity/pcoa/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "pcoa",
                      dist.metric,host,
                      comparison,agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=pcoa.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = image.format)
}



# PCA ####
# Convert df into wide format
ps.q.df.wide.pca<-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()


# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide.pca)<-ps.q.df.wide.pca$Sample
ps.q.df.wide.pca<-ps.q.df.wide.pca[,-1] # remove sample names
ps.q.df.wide.pca<-as.matrix(ps.q.df.wide.pca) # convert to matrix
# The first crucial step is to transform the compositional data:
# Use a log-ratio transformation, typically the centered log-ratio (clr)
# or isometric log-ratio (ilr) transformation.
# The clr transformation is often preferred for interpretability, 
# while ilr is used for statistical properties.
# We will use robust CLR (RCLR) with natural log
ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,method="rclr",logbase = exp)

#### >if you want to exclude specific samples
if(host=="NMR"){
  ps.q.pruned<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("M40")),]
}else if(host=="mice"){
  ps.q.pruned<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("mf_1","MSM343")),]
}
# remove ASVs that are zero because their samples are removed
ps.q.pruned<-ps.q.pruned[,which(colSums(ps.q.pruned)!=0)]
ps.q.df.wide.pca.tfm<-decostand(ps.q.pruned,method="rclr",logbase = exp)
####<

# calculate principal components
pca.q<-prcomp(ps.q.df.wide.pca.tfm)
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
  left_join(custom.md)


if(comparison=="age"){
  pca.plot<-pca.plot%>%
    ggplot(aes(x=PC1,y=PC2,color=age,fill=age))+
    scale_fill_viridis_c(option = "D")+
    labs(fill="Age (years)")
}else if (comparison=="sex"){
  pca.plot<-pca.plot%>%
    ggplot(aes(x=PC1,y=PC2,color=sex,fill=sex))+
    scale_color_manual(breaks = custom.levels,
                       labels=unname(pretty.level.names),
                       values = custom.colors)+
    scale_fill_manual(name=NULL, breaks = custom.levels,
                      labels=custom.levels,
                      values = custom.colors)+
    guides(fill="none")+
    stat_ellipse(geom = "polygon",
                 level = 0.8,
                 alpha=0.2,
                 show.legend = FALSE)
}else if(comparison=="strain"){
  pca.plot<-pca.plot%>%
    ggplot(aes(x=PC1,y=PC2,color=class,fill=class))+
    scale_color_manual(breaks = custom.levels,
                       labels=unname(pretty.level.names),
                       values = custom.colors)+
    scale_fill_manual(name=NULL, breaks = custom.levels,
                      labels=custom.levels,
                      values = custom.colors)+
    guides(fill="none")+
    stat_ellipse(geom = "polygon",
                 level = 0.8,
                 alpha=0.2,
                 show.legend = FALSE)
}

pca.plot<-pca.plot  +
  geom_point(size=3,pch=21)+
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"),
       color=gg.labs.name)+
  theme_classic()+
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  ggtitle(paste("PCA between different",as.character(host.class[host]),gg.title.groups, "\n",
                gg.title.taxon))+
  scatter.plot.theme  # theme for scatter plots (PCA, PCoA, nMDS)

library(ggrepel)
pca.labeled<-pca.plot+
  ggrepel::geom_text_repel(aes(label=Sample),
                           show.legend = FALSE,size=7) # add labels to samples

if(comparison=="age"){
  pca.labeled<-pca.plot+
    ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),max.overlaps = Inf,
                             show.legend = FALSE,size=7) # add labels to samples
}

for(image.format in image.formats){
  ggsave(paste0("./images/diversity/pca/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "pca",
                      host,
                      comparison,agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=pca.plot,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = image.format)
  ggsave(paste0("./images/diversity/pca/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "pca-labeled",
                      host,
                      comparison,agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=pca.labeled,
         width = 4500,height = 3000,
         units = "px",dpi=300,device = image.format)
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


foo<-all.div%>%
  select(-tot,-sex)%>%
  pivot_wider(names_from = metric)%>%
  left_join(ps.q.agg)

div.otus<-foo%>%
  filter(agegroup=="[0,10)")%>%
  select(Sample,sobs,OTU,Abundance)%>%
  pivot_wider(names_from = OTU,
              values_from = Abundance,
              values_fill = 0)%>%
  distinct(Sample,.keep_all = T)
cov.matr<-cov(div.otus[,-1])
cov.matr[-1,1]%>%
  as.data.frame()%>%
  rename("sobs"=".")%>%
  rownames_to_column("OTU")%>%
  as_tibble()%>%
  arrange(sobs)

foo%>%
  filter(OTU%in%c("5a818007a6452d5db9223ef1388e451c",
                  "2a9b76ddc8fc1e5f0e10f1fd40de1c79",
                  "c513b4ef037af47125244161fb1eab50",
                  "72aa784b6a3f05f45910fe33902caa81",
                  "f522ba24c2cefcb9eaecc3bd442fc926",
                  "2d4fb4e16e2ec5372a78961fcc33e48c",
                  "65441a076e1d645106982a3e2423e0b3",
                  "cbb7a2fe6431a3d874d38a8b6a104854"))%>%
  ggplot(aes(x=Abundance,y=sobs))+
  geom_point()+
  facet_wrap(~OTU)

foo%>%
  filter(OTU%in%c("5a818007a6452d5db9223ef1388e451c",
                  "2a9b76ddc8fc1e5f0e10f1fd40de1c79",
                  "c513b4ef037af47125244161fb1eab50",
                  "72aa784b6a3f05f45910fe33902caa81",
                  "f522ba24c2cefcb9eaecc3bd442fc926",
                  "2d4fb4e16e2ec5372a78961fcc33e48c",
                  "65441a076e1d645106982a3e2423e0b3",
                  "cbb7a2fe6431a3d874d38a8b6a104854"))%>%
  ungroup%>%
  select(Genus)%>%
  distinct(Genus)
