library(tidyverse)
library(phyloseq)
library(indicspecies)

truncationlvl<-"234"
agglom.rank<-"Genus"

read.end.type<-"single"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))


pretty.facet.labels<-
  c(
    # "NMR" = "*Heterocephalus glaber*", # better labels for facets
    "B6mouse" = "B6 mouse",
    "MSMmouse" = "MSM/Ms mouse",
    "FVBNmouse" = "FVB/N mouse"
    # ,
    # "DMR" = "*Fukomys Damarensis*",
    # "hare" = "*Lepus europaeus*",
    # "rabbit" = "*Oryctolagus cuniculus*",
    # "spalax" = "*Nannospalax leucodon*",
    # "pvo" = "*Pteromys volans orii*"
  )

custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]


ps.q.df <-ps.q.agg.abs%>%
  filter(class%in%custom.levels)%>%
  # filter(Taxon.bp!="Remainder (Mean abundance < 1%)")%>%
  select(Sample,OTU,Abundance,class,Taxon)

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
abund<-ps.q.df.wide[,-c(1,2)]  
sample.info<-ps.q.df.wide$Sample
class.info<-ps.q.df.wide$class

# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  
ps.q.df.wide<-as.matrix(ps.q.df.wide)

# Run indicator species command
set.seed(1)
inv = multipatt(abund, class.info, func = "r.g", control = how(nperm=9999))

summary(inv)

ps.q.agg.abs%>%
  filter(Taxon=="Parasutterella (Sutterellaceae)",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Parasutterella (Sutterellaceae)")

# cors<-cor(abund,method="spearman")

