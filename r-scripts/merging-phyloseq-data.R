# merging phyloseq data
# load biagi data
authorname<-"biagi"
agglom.rank<-"Genus"
truncationlvl<-"0"
read.end.type<-"single"
load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))
ps.b<-ps.q
ps.b.agg<-ps.q.agg
ps.b.total<-ps.q.total
ps.b.1pc<-ps.q.1pc
custom.md.b<-custom.md

# load pooled data
authorname<-"pooled"
agglom.rank<-"Genus"
truncationlvl<-"234"
read.end.type<-"single"

load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))
# merge metadata
custom.md<-rbind(custom.md,custom.md.b)

# merge all data
ps.q.agg<-rbind(ps.q.agg,ps.b.agg)
ps.q.total<-rbind(ps.q.total,ps.b.total)
ps.q.1pc<-rbind(ps.q.1pc,ps.b.1pc)

rm(ps.b)
rm(ps.b.agg)
rm(ps.b.1pc)
rm(ps.b.total)
rm(custom.md.b)
authorname<-"merged"
directory<-gsub("pooled","merged",directory)
qiimedir<-gsub("pooled","merged",qiimedir)
save.image(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                      truncationlvl,agglom.rank,
                                      "phyloseq-workspace.RData",sep = "-")))
