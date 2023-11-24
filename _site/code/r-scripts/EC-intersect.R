nmr.sim<-c("EC:1.17.1","EC:1.3.99","EC:1.97.1","EC:3.1.21","EC:3.2.1",
           "EC:3.5.99","EC:3.6.1")
nmr.perm.full<-c("EC:3.1.4","EC:1.8.4","EC:1.1.1","EC:3.4.22","EC:1.3.5",
                 "EC:1.20.4","EC:1.2.7","EC:1.6.5","EC:1.4.1", "EC:3.4.21",
                 "EC:3.2.2","EC:3.4.13","EC:1.17.4","EC:3.3.1","EC:3.2.1",
                 "EC:3.5.99")

nmr.perm.fast<-c("EC:3.1.4","EC:1.8.4","EC:1.5.1","EC:3.1.26")

intersect(nmr.sim,nmr.perm.full)
intersect(nmr.sim,nmr.perm.fast)
intersect(nmr.perm.full,nmr.perm.fast)


spfmouse.sim<-c("EC:1.15.1","EC:1.16.3","EC:1.17.1","EC:1.2.7","EC:1.3.98",
                "EC:1.3.99","EC:1.4.1","EC:1.97.1","EC:3.1.1","EC:3.1.21",
                "EC:3.1.22","EC:3.11.1","EC:3.2.1","EC:3.2.2","EC:3.4.11",
                "EC:3.4.19","EC:3.4.21","EC:3.5.1","EC:3.5.2","EC:3.5.99",
                "EC:3.6.1")


spfmouse.perm.full<-c("EC:1.17.4","EC:1.3.99","EC:1.2.1","EC:1.3.5","EC:1.4.1",
                      "EC:3.2.2","EC:1.8.1","EC:3.4.22","EC:1.97.1","EC:3.3.1",
                      "EC:1.12.99","EC:3.1.5","EC:3.1.21","EC:1.8.4","EC:1.18.1",
                      "EC:3.2.1","EC:3.4.13","EC:1.12.1","EC:3.1.13",
                      "EC:1.7.2","EC:3.4.19","EC:3.4.21")
intersect(spfmouse.sim,spfmouse.perm.full)

intersect(nmr.perm.full,spfmouse.perm.full)

common.enz<-intersect(spfmouse.perm.full,nmr.perm.full)
nmr.enz<-setdiff(nmr.perm.full,spfmouse.perm.full)
spfmouse.enz<-setdiff(spfmouse.perm.full,nmr.perm.full)

knitr::kable(common.enz,format = "simple")
knitr::kable(nmr.enz,format = "simple")
knitr::kable(spfmouse.enz,format = "simple")

filtered.enzymes.long<-filtered.enzymes.long%>%
  inner_join(metadata,by=c("Sample"="sample.id"))

filtered.enzymes.long%>%
  # filter(class=="NMR")%>%
  filter(predicted.function%in%names(signif.enz))%>%
  group_by(predicted.function,class)%>%
  summarise(n=n())%>%View

filtered.enzymes.long%>%
  group_by(predicted.function)%>%
  # filter(predicted.function%in%nmr.perm.full)%>%
  mutate(class.subset=ifelse(predicted.function%in%common.enz, "common",
                       ifelse(predicted.function%in%nmr.enz, "unique_for_NMR",
                              ifelse(predicted.function%in%spfmouse.enz, "unique_for_SPFmouse",
                                     "other"))))%>%
  ggplot(aes(x=class.subset,y=Abundance,
         fill=factor(predicted.function)))+
  geom_bar(stat = "identity")
