dist.metric<-"bray"
dist.metric<-"jaccard"
dist.metric<-"aitchison"
dist.metric<-"canberra"

ps.sampledata<-ps.q.agg.abs%>%
  select(c("Sample","class","origin","relation",
           "sex","caste","birthday"))%>%
  distinct() # metadata

ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)
all.wide<-ps.q.df%>%
  select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs and rownames are sample IDs
rownames(all.wide)<-all.wide$Sample
all.wide<-all.wide[,-c(1,2)]  
all.wide<-as.matrix(all.wide)

a.size<-rowSums(all.wide)
a.shape<-all.wide/a.size
a.shape
if(dist.metric=="jaccard"){
  dist<-avgdist(all.wide,
                dmethod="jaccard",
                binary=TRUE,
                sample=1000,
                iterations = 1000)
}else if(dist.metric=="bray"|dist.metric=="canberra"){
  dist<-avgdist(all.wide,
                dmethod=dist.metric,
                sample=1000,
                iterations = 1000)
}else if(dist.metric=="aitchison"){
  # dist.tfm<-decostand(x = all.wide,"clr",pseudocount=1)
  dist<-avgdist(all.wide,
                dmethod="aitchison",
                sample=1000,
                iterations = 1000,
                pseudocount=1)
}
dist<-vegdist(all.wide,method = dist.metric)
dist<-vegdist(all.wide,method = "chisq")

range(a.size)


a.shape<-as.data.frame(a.shape)
# adonis.AD<-adonis2(dist~.+a.size,data=a.shape,permutations = 0)
# nterms<-length(adonis.AD$aov.tab$R2)

# A =size not shape, B=size and shape, C=shape not size, D=residual

adonis.AD<-adonis(dist~.+a.size,data=a.shape,permutations = 0)
nterms<-length(adonis.AD$aov.tab$R2)
AD<-adonis.AD$aov.tab$R2[c(nterms-2,nterms-1)]
A<-AD[1]
D<-AD[2]

# C is shape when size is partialled out and without the residual
adonis.C<-adonis(dist~a.size+.,data = a.shape,permutations = 0)
nterms<-length(adonis.C$aov.tab$R2)
C<-1-sum(adonis.C$aov.tab$R2[c(1,nterms-1)])

B<-1-A-C-D
c(A,B,C,D)*100
barplot(c(A,B,C,D)*100,names.arg = c(A,B,C,D)*100)

# when variables>samples, there's no pure size and no residual

# To take both size and shape into account jointly
# Euclidean on log tfmed? log(1+x)
# chisq on raw?




a<-c(7,6,4,0)
b<-a*10
df.aa<-t(data.frame(a,a))
df.ab<-t(data.frame(a,b))

vegdist(df.ab) # bray 0
vegdist(df.ab) # 0.8181818

vegdist(df.aa,method = "jaccard",binary = T) # jaccard 0
vegdist(df.ab,method = "jaccard",binary = T)

vegdist(df.aa,method = "canberra") # canberra 0
vegdist(df.ab,method = "canberra") # 0.8181818

vegdist(df.aa,method = "aitchison",pseudocount=1) # aitchison 0
vegdist(df.ab,method = "aitchison",pseudocount=1) # 0

vegdist(df.aa,method = "euclidian") # euclidian 0
vegdist(df.ab,method = "euclidian") # 90.44888

decostand(df.ab,method="clr",pseudocount=1)

library(DESeq2)
vst.matrix<-varianceStabilizingTransformation(t(all.wide),fitType = "local")%>%
  t()

mat<-t(as.matrix(df.ab)+1)
varianceStabilizingTransformation(mat,fitType = "local")

vst.dist<-vegdist(vst.matrix,method="euclidian")%>%
  as.matrix()%>%
  as_tibble(rownames="Sample")%>%
  pivot_longer(-Sample)%>%
  filter(name<Sample)



