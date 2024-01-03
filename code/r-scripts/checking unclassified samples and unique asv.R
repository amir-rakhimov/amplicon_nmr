ps.q.agg[grep("Unclassified|Uncultured",ps.q.agg$Taxon),]%>%
  filter(Abundance>0)%>%
  group_by(Sample)%>%
  summarise(n=sum(Abundance))%>%
  arrange(-n)%>%
  inner_join(.,unique(ps.q.agg[,c("class","Sample")]))%>%
  filter(n>30)%>%
  group_by(class)%>%
  summarise(nh=n())
View()


ps.q.agg%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(OTU)%>%
  summarise(OTU_count=n())%>%
  arrange(-OTU_count)

ps.q.agg%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(Taxon)%>%
  summarise(Taxon_count=n())%>%
  arrange(-Taxon_count)
