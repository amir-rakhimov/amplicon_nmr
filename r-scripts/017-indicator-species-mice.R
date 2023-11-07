library(tidyverse)
library(phyloseq)
library(indicspecies)
library(vegan)

truncationlvl<-"234"
agglom.rank<-"OTU"

read.end.type<-"single"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

rare.status<-"nonrarefied"
filter.status<-"nonfiltered"

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
  filter(Taxon=="Unclassified (Tannerellaceae)",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Unclassified (Tannerellaceae)")

# cors<-cor(abund,method="spearman")

# Save output to a file
sink(paste0("./output/multipatt-",rare.status,"-",filter.status,"-",agglom.rank,
            "-",paste(custom.levels,collapse = "-"),".txt"))
cat(paste(custom.levels,collapse = " "),"\n")
summary(inv)
sink()

# Rarefied
min.n_seqs.all<-ps.q.agg.abs%>%
  filter(class %in% custom.levels)%>%
  select(Sample,OTU,Abundance)%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# rarefied asv table with vegan
set.seed(1)
if(agglom.rank!="OTU"){
  ps.q.df.rare<-rrarefy(ps.q.df.wide,sample=min.n_seqs.all)%>%
    as_tibble(rownames="Sample")%>%
    pivot_longer(-Sample)%>%
    as.data.frame()%>%
    inner_join(unique(ps.q.agg.abs[,c("Sample","class")]),
               by="Sample")%>%
    rename(Taxon=name,Abundance=value)
}else{
  ps.q.df.rare<-rrarefy(ps.q.df.wide,sample=min.n_seqs.all)%>%
    as_tibble(rownames="Sample")%>%
    pivot_longer(-Sample)%>%
    as.data.frame()%>%
    inner_join(unique(ps.q.agg.abs[,c("Sample","class")]),
               by="Sample")%>%
    rename(OTU=name,Abundance=value)
}

if(agglom.rank!="OTU"){
  ps.q.df.rare.wide<-ps.q.df.rare%>%
    pivot_wider(names_from = "Taxon",
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()
}else{
  ps.q.df.rare.wide<-ps.q.df.rare%>%
    pivot_wider(names_from = "OTU", 
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()
}
abund.rare<-ps.q.df.rare.wide[,-c(1,2)]  
sample.info.rare<-ps.q.df.rare.wide$Sample
class.info.rare<-ps.q.df.rare.wide$class

save.image(paste0("./rdafiles/indicspecies-",paste(custom.levels,collapse = "-"),
                  "-",agglom.rank,"-workspace.RData"))
# Run indicator species command
set.seed(1)
inv.rare = multipatt(abund.rare, sample.info.rare, func = "r.g", control = how(nperm=9999))

summary(inv.rare)


