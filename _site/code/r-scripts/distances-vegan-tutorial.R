library(tidyverse)
library(vegan)
days_wanted<-c(0:9,141:150)
shared<-read_tsv("mice.shared") %>%
  select(Group,starts_with("Otu"))%>%
  mutate(day=str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day)%>%
  pivot_longer(-Group) %>% 
  group_by(Group)%>%
  mutate(total=sum(value)) %>%
  filter(total>1800) %>%
  group_by(name) %>%
  mutate(total=sum(value)) %>%
  filter(total!=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Group) %>% 
  as.data.frame()

rownames(shared)<-shared$Group
shared<-shared[,-1]
shared<-as.matrix(shared)

# distances with vegan ####
dist<-vegdist(shared, method="bray")
set.seed(19760620)
nmds<-metaMDS(dist)

scores(nmds)%>%
  as.tibble(rownames="Group")%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point()+
  theme_bw()

# Using community matrix directly in metaMDS
nmds<-metaMDS(shared,autotransform = FALSE)
# doesn't work tbh


# rarefaction ####
set.seed(1)
dist<-avgdist(shared,dmethod="bray",sample=1800)
set.seed(1)
nmds<-metaMDS(dist)
scores(nmds)%>%
  as.tibble(rownames="Group")%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point()+
  theme_bw()



# envfit ####
data(varespec, varechem)
# varespec has 44 species for 24 sites
# varechem has soil characteristics for the same 24 sites
library(MASS)
ord <- metaMDS(varespec)
(fit <- envfit(ord, varechem, perm = 999))
scores(fit, "vectors")
plot(ord)
plot(fit)
plot(fit, p.max = 0.05, col = "red")

## Adding fitted arrows to CCA. We use "lc" scores, and hope
## that arrows are scaled similarly in cca and envfit plots
ord <- cca(varespec ~ Al + P + K, varechem)
plot(ord, type="p")
fit <- envfit(ord, varechem, perm = 999, display = "lc")
plot(fit, p.max = 0.05, col = "red")

## 'scaling' must be set similarly in envfit and in ordination plot
plot(ord, type = "p", scaling = "sites")
fit <- envfit(ord, varechem, perm = 0, display = "lc", scaling = "sites")
plot(fit, col = "red")




#####
myTable = matrix(nrow=20,ncol=10)
for(i in 1:10) {myTable[i,] = sample(1:10,10)}
for(i in 11:20) {myTable[i,] = sample(15:25,10)}

myMetadata = data.frame(FirstColumn = character(), SecondColumn = character(), stringsAsFactors=FALSE)
myMetadata[1:10,1] = 'G1'
myMetadata[11:20,1] = 'G2'
myMetadata[seq(1,20,2),2] = 'G3'
myMetadata[seq(2,20,2),2] = 'G4'

envfit(myTable ~ myMetadata$FirstColumn, data = myMetadata, perm=1000)

envfit(myTable ~ myMetadata$SecondColumn, data = myMetadata, perm=1000)


sites=c("Site A","Site B","Site C","Site D","Site E","Site F","Site 
G","Site H","Site I","Site J","Site K","Site L","Site M","Site N","Site O","Site P","Site Q","Site R","Site S","Site T","Site U")
american.elm=c(41.91,10.11,2.62,5.31,7.51,9.72,17.44,9.06,19.83,30.81,62.6,21.29,20.7,28.68,27.69,34.89,35.65,3.87,12.68,1.58,2.97)
white.birch=c(7.07,15.89,26.77,15.61,14.59,6.33,2.23,11.66,21.49,20.15,7.61,23.29,0,0,0,0,0,0,0,56.09,42.34)
red.oak=c(0,0,0,0,0,0,0,0,0,0,0,6.02,0,0,0,0,0,0,0,0,0.05)
populus.grand=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.11,0)
beech=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.36,5.45)
sugar.maple=c(0.49,2.64,3.35,4.6,3.37,2,1.32,4.21,4.13,3.61,0.34,1.2,0,0,0,0,0,0,0,2.19,0.09)

mytable <- data.frame(sites,american.elm,red.oak,populus.grand,beech,sugar.maple)
mytable<-mytable[,2:ncol(mytable)]

mytable.NMDS=metaMDS(mytable, distance = "jaccard", k = 4, trymax = 2000, autotransform=FALSE)

plot.mytable<-data.frame(mytable.NMDS$points)
plot.mytable
par(mar=c(3,3,2,5) ,mgp=c(1.8,0.75,0))
plot(plot.mytable$MDS1, plot.mytable$MDS2, pch=16, cex=1, col="black",
     xlab="NMDS1", ylab="NMDS2", cex.lab=1, cex.axis=1, main="", bty="L",
     mai=c(0,0,2,10), xlim=c(-1.5,1.3), ylim=c(-0.9,1))

fit <- envfit(mytable.NMDS, mytable, choices=c(1,2,3))
fit.plot = plot(fit, cex=1.3, col="red", xlim=c(-1.5,1.3), ylim=c(-1.2,1.2),
                xlab="NMDS1",ylab="NMDS2")
fit

Trees=c("american.elm","red.oak","populus.grand","beech","sugar.maple")
Tree.NMDS1=c(-0.76538,-0.1533,0.36065,0.25411,0.49583)
Tree.NMDS2=c(-0.27961,0.06605,-0.51345,-0.79497,0.84299)
Tree.NMDS.scores=data.frame(Trees,Tree.NMDS1,Tree.NMDS2)
# Overlay the NMDS score on the plot
# These values are way far from what they should be
points(Tree.NMDS.scores$Tree.NMDS1,Tree.NMDS.scores$Tree.NMDS2,
       col="red", pch=16)

fit$vectors$r
ordiArrowMul(fit)

scrs<-scores(fit, "vectors", choices= 1:2)
scrs
scrs * ordiArrowMul(fit)

plot(mytable.NMDS, display = "sites", type = "n")
points(mytable.NMDS, display = "sites", pch = 19, col = "black")
plot(fit, col = 'red')

## add the locations of arrow heads as blue points to see if the correspond
points(scrs * ordiArrowMul(fit), col = "blue")

