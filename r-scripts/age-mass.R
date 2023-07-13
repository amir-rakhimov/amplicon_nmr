library(ggplot2)
library(ggrepel)
workdir<-"./mega/data/"
lsp<-read.table(paste0(workdir,"lifespans and body traits.csv"),
                header = T,sep ='\t')
str(lsp)
colnames(lsp)<-c("taxon","lifespan","body.mass.kg",
                 "avg.brain.size.kg","Remarks")
dat<-lsp
dat<-dat[!grepl("Not yet established",dat$lifespan),]
dat<-dat[!is.na(dat$lifespan),]
dat<-dat[!is.na(dat$body.mass.kg),]
dat<-dat[!grepl("Homo sapiens",dat$taxon),]
dat$lifespan<-as.numeric(dat$lifespan)
dat$loglif<-log(dat$lifespan)
dat$logmas<-log(dat$body.mass.kg)

agemasplot<-ggplot(dat,aes(x=logmas,y=loglif))+
  geom_point(col=2)+
  geom_smooth(method = "lm",se=FALSE)+ # fit linear model
  geom_text_repel(aes(x=logmas,y=loglif,
                      label=taxon))+
  theme_bw()+
  xlab("log Body mass (kg)")+
  ylab("log Longevity (years)")+
  ggtitle("The relationship between body mass and longevity for rodents")
  

ggsave(filename=paste0("./images/lineplots/",
                       Sys.Date(),"-age-mass.png"),agemasplot,
      width = 2000,height = 2000, units="px",
      dpi = 300)
