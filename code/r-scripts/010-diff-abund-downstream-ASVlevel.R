library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)
truncationlvl<-"234"
agglom.rank<-"OTU"
authorname<-"pooled"

rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"

host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
comparison<-"sex"
comparison<-"strain"

host.labels<-
  c("NMR" = "*Heterocephalus glaber*")

host.labels<-
  c("B6mouse" = "B6 mouse",
    "MSMmouse" = "MSM/Ms mouse",
    "FVBNmouse" = "FVB/N mouse")

ref.level<-"MSMmouse"
load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))
maaslin.signif.features<-list()
aldex.signif.features<-list()
ancombc.signif.features<-list()
for (ref.level in custom.levels){
  lvl.maaslin.signif.features<-
    read.table(file.path("./rtables/pooled",
                         paste("maaslin",rare.status,
                               filter.status,host,agglom.rank,
                               comparison,truncationlvl,ref.level,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  lvl.aldex.signif.features<-
    read.table(file.path("./rtables/pooled",
                         paste("aldex2",rare.status,
                               filter.status,host,agglom.rank,
                               comparison,truncationlvl,ref.level,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  lvl.ancombc.signif.features<-
    read.table(file.path("./rtables/pooled",
                         paste("ancombc",rare.status,
                               filter.status,host,agglom.rank,
                               comparison,truncationlvl,ref.level,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  lvl.ancombc.signif.features<-lvl.ancombc.signif.features%>%
    dplyr::select(-X.Intercept.)%>%
    pivot_longer(!taxon_id,names_to = "class",values_to = "coef")
  
  maaslin.signif.features[[ref.level]]<-lvl.maaslin.signif.features
  aldex.signif.features[[ref.level]]<-lvl.aldex.signif.features
  ancombc.signif.features[[ref.level]]<-lvl.ancombc.signif.features
}
rm(lvl.maaslin.signif.features)
rm(lvl.aldex.signif.features)
rm(lvl.ancombc.signif.features)

maaslin.signif.features<-bind_rows(maaslin.signif.features,.id = comparison)
aldex.signif.features<-bind_rows(aldex.signif.features,.id = comparison)
ancombc.signif.features<-bind_rows(ancombc.signif.features,.id = comparison)

maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by(strain,feature)%>%
  filter(n()==length(custom.levels)-1)%>% 
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))

aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)

ancombc.signif.decreased<-ancombc.signif.features%>%
  as_tibble()%>%
  group_by(strain,taxon_id)%>%
  summarize(count=sum(coef<0))%>% # coef should be negative in all other groups
  # so we count how many coefs are negative per taxon_id
  filter(count==length(custom.levels)-1)%>%
  left_join(ancombc.signif.features,by=c("strain"="strain", "taxon_id"="taxon_id"))

Reduce(intersect,list(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="B6mouse"],
                      aldex.neg.effect$OTU[aldex.neg.effect$strain=="B6mouse"],
                      ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="B6mouse"]))

# ALDEX shares no otu with other two datasets simultaneously
# but separately, yes
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="B6mouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="B6mouse"])
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="B6mouse"],
          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="B6mouse"])
intersect(ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="B6mouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="B6mouse"])

intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="MSMmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="MSMmouse"])
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="MSMmouse"],
          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="MSMmouse"])
intersect(ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="MSMmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="MSMmouse"])

intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="FVBNmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="FVBNmouse"])
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="FVBNmouse"],
          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="FVBNmouse"])
intersect(ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="FVBNmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="FVBNmouse"])

signif.all.groups<-list()
for (ref.level in custom.levels){
  signif.lvl<-
    Reduce(intersect,list(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain==ref.level],
                          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain==ref.level]))
  # filter by mean relative abundance: >=0.5%
  signif.lvl.filtered<-ps.q.agg%>%
    group_by(OTU,class)%>%
    filter(OTU %in% signif.lvl)%>%
    dplyr::select(OTU, MeanRelativeAbundance)%>%
    arrange(-MeanRelativeAbundance)%>%
    distinct(OTU,.keep_all = T)%>%
    ungroup()%>%
    pull(OTU)%>%
    unique()
  signif.all.groups[[ref.level]]<-signif.lvl.filtered
  write.table(signif.lvl.filtered,
              file=file.path("./rtables/pooled",
                             paste("significant-features",rare.status,
                                   filter.status,host,agglom.rank,
                                   comparison,truncationlvl,
                                   ref.level,"signif.tsv",sep="-")),
              row.names = F,sep = "\t",col.names = F)
}
rm(signif.lvl)

save(signif.all.groups,
     file=paste0("./rdafiles/",
               paste("significant-features-all-groups",
                     rare.status,filter.status,host,agglom.rank,
                     comparison,truncationlvl,
                     sep="-"),".RData"))
