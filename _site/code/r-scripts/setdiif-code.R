ps.q.agg.rel<-ps.q %>%
  #subset_samples(class=="NMR") %>% # choose only naked mole-rat
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  transform_sample_counts(function(x)100* x / sum(x)) %>% # transform into percentages
  psmelt()  # 

nmr.foo<-ps.q.agg.rel%>%filter(class=="NMR",Abundance!=0)
length(table(nmr.foo$OTU))

spf.foo<-ps.q.agg.rel%>%filter(class=="SPFmouse",Abundance!=0)
length(table(spf.foo$OTU))

FukomysDamarensis.foo<-ps.q.agg.rel%>%filter(class=="FukomysDamarensis",Abundance!=0)
length(table(FukomysDamarensis.foo$OTU))

pal.foo<-ps.q.agg.rel%>%filter(class=="pal",Abundance!=0)
length(table(pal.foo$OTU))

setdiff(names(table(nmr.foo$Taxon.bp)),names(table(spf.foo$Taxon.bp)))


