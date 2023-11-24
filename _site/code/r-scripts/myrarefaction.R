ps.sampledata<-ps.q.agg.abs%>%
  select(c("Sample","class","origin","relation",
           "sex","caste","birthday"))%>%
  distinct() # metadata

ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)
all.wide<-ps.q.df%>%
  select(-OTU,-class)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

all.wide<-all.wide%>%pivot_longer(-Sample)

# null model: repeat each row by the value in the 'value' column
# e.g otu1 will be shown 115 times
# then, randomise the order of otus and recount everything
rand<-all.wide%>%
  uncount(value)%>% # to repeat each row and get 1.8 mln rows
  mutate(rand_name=sample(name)) %>% # randomise without replacement
  select(-name)%>%
  dplyr::count(Sample, rand_name)

# sum of reads for each sample
rand_group_count<-rand%>%
  group_by(Sample)%>%
  summarise(n=sum(n))
original_group_count<-all.wide%>%
  group_by(Sample)%>%
  summarise(n=sum(value))

inner_join(original_group_count,rand_group_count,by="Sample")

# sum of reads for each otu
original_otu_count<-all.wide%>%
  group_by(name)%>%
  summarise(n=sum(value))
rand_otu_count<-rand%>%
  group_by(rand_name)%>%
  summarise(n=sum(n))
inner_join(original_otu_count,rand_otu_count,by=c("name"="rand_name"))


# each sample is a statistical sample of the overall community distro

# we preserve the total number of sequences in each otu and 
# total num of sequences in each sample

rand_group_count%>%
  ggplot(aes(x=n))+
  geom_histogram()

# turn rand into matrix
rand_df<-rand%>%
  pivot_wider(names_from = "rand_name",
              values_from = "n",
              values_fill = 0)%>%
  as.data.frame()
rownames(rand_df)<-rand_df$Sample
rand_df<-rand_df[,-1]

rand_matrix<-as.matrix(rand_df)

# no rarefaction distance matrix
norare_dist_matrix<-vegdist(rand_matrix,method="bray")
# rarefaction distance matrix
set.seed(1)
rare_dist_matrix<-avgdist(rand_matrix,
                          dmethod="bray",
                          sample=min.n_seqs.all) # size is the smallest library size

# compare distance calculation approaches
norare_dist_tibble<-norare_dist_matrix%>%
  as.matrix()%>%
  as_tibble(rownames="Sample")%>%
  pivot_longer(-Sample)%>%
  filter(name<Sample) # only unique comparisons (lower triangle)


rare_dist_tibble<-rare_dist_matrix%>%
  as.matrix()%>%
  as_tibble(rownames="Sample")%>%
  pivot_longer(-Sample)%>%
  filter(name<Sample) # only unique comparisons (lower triangle)

# range of rarefied data is small, while for non-rarefied it's huge: 0-1
# rarefied distances are averaging about the same
comparison<-inner_join(norare_dist_tibble,rare_dist_tibble,by=c("Sample","name"))%>%
  select(Sample,name,norare=value.x,rare=value.y)%>%
  inner_join(.,rand_group_count,by=c("Sample"="Sample"))%>% # add total counts of left sample
  inner_join(.,rand_group_count,by=c("name"="Sample"))%>% # add total counts of right sample (name)
  mutate(n_diff=abs(n.x-n.y))%>%
  select(-n.x,-n.y)

comparison%>%
  ggplot(aes(x=norare,y=rare,color=n_diff))+
  geom_point(size=0.25,alpha=0.25)+
  geom_smooth()

# make a faceted plot where x axis is the number of sequences that's 
# different between two compared samples, and y axis is distances
comparison%>%
  # split the df: for each comparison, there's a row for dist in rare
  # and one row for dist in norare
  pivot_longer(cols=c("norare","rare"),names_to="type",values_to = "dist")%>%
  ggplot(aes(x=n_diff,y=dist))+
  geom_point(size=0.25,alpha=0.25)+
  facet_wrap(~type,nrow = 2)

# rarefaction removes the variability due to sequencing depth


dist<-vegdist(all.wide,method = "bray")
set.seed(1)
nmds<-metaMDS(dist)
# nmds<-metaMDS(shared,autotransform = FALSE)

scores(nmds)%>%
  as_tibble(rownames="Sample")%>% # convert rownames into "Sample" column
  inner_join(.,ps.sampledata, by="Sample")%>%
  ggplot(aes(x=NMDS1,y=NMDS2,color=class))+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)

set.seed(1)
dist<-avgdist(all.wide,dmethod = "bray",sample=1000,iterations = 1000)
set.seed(1)
nmds<-metaMDS(dist)
scores(nmds)%>%
  as_tibble(rownames="Sample")%>% # convert rownames into "Sample" column
  inner_join(.,ps.sampledata, by="Sample")%>%
  ggplot(aes(x=NMDS1,y=NMDS2,color=class))+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  geom_point()+
  scale_color_manual(breaks = custom.levels,
                     labels=custom.color.labels,
                     values = custom.colors)
