library(tidyverse)
directory<-"./data/yasuda-data/"
workdir<-"/home/amir/yasuda/data/"
filenames<-read.table(paste0(directory,
                             "filenames-yasuda-v1.csv"),
                      header = T,
                      sep = "\t")
srr.list<-read.table(paste0(directory,
                            "SraRunTable.txt"),
                     header = T,
                     sep ="\t")
head(filenames)
colnames(srr.list)<-c("sample")

# fwd.reads<-filenames$V1[grep("pass_1",filenames$V1)]
# rev.reads<-filenames$V1[grep("pass_2",filenames$V1)]

fwd.reads<-paste0(filenames$Run,"_1.fastq.gz")
rev.reads<-paste0(filenames$Run,"_2.fastq.gz")
filenames$forward<-paste0(workdir,fwd.reads)
filenames$reverse<-paste0(workdir,rev.reads) 

colnames(filenames)<-c("sample-id",
                       "host",
                       "host-body-site",
                       "host-subject-id",
                       "sex",
                       "forward-absolute-filepath",
                      "reverse-absolute-filepath")

write.table(filenames,
            file = paste0(directory,"filenames-anders.tsv"),
            sep = "\t",row.names = F, quote = F)





####
directory<-"./data/yasuda-data/"
workdir<-"/home/amir/yasuda/data/"
metadata<-read.table(paste0(directory,
                             "filenames-yasuda-v1.csv"),
                      header = T,
                      sep = "\t")
filenames<-list.files(paste0(directory,"fastq"))
fwd.reads<-filenames[grep("R1_001.fastq.gz",filenames)]
rev.reads<-filenames[grep("R2_001.fastq.gz",filenames)]
metadata$forward<-paste0(workdir,fwd.reads)
metadata$reverse<-paste0(workdir,rev.reads)
View(metadata[,c("sample.id","forward","reverse")])
write.table(metadata,
            file = paste0(directory,"filenames-yasuda-v2.tsv"),
            sep = "\t",row.names = F, quote = T)


####
directory<-"./data/shanmuganandam-data/metadata/"
workdir<-"/home/amir/shanmuganandam/data/final/"
metadata<-read.table(paste0(directory,
                            "filenames-shanmuganandam.tsv"),
                     header = T,
                     sep = "\t")

fwd<-metadata[seq(from=1, to=nrow(metadata),by=2),]
rev<-metadata[seq(from=2, to=nrow(metadata),by=2),]

table(fwd$sex==rev$sex)
md<-cbind(fwd,rev$filename)
md$filename<-paste0(workdir,md$filename)
md$`rev$filename`<-paste0(workdir,md$`rev$filename`)
md<-md%>%
  arrange(strain)
colnames(md)<-c("sample-id",
                "animal",
                "strain",
                "sex",
                "forward-absolute-filepath",
                "reverse-absolute-filepath")
write.table(md,
            file = paste0(directory,"filenames-shanmuganandam-local.tsv"),
            sep = "\t",row.names = F, quote = F)


####
directory<-"./data/liu-data/metadata/"
workdir<-"/home/amir/liu/data/final/"
metadata<-read.table(paste0(directory,
                            "filenames-liu.tsv"),
                     header = T,
                     sep = "\t")

fwd<-metadata[seq(from=1, to=nrow(metadata),by=2),]
rev<-metadata[seq(from=2, to=nrow(metadata),by=2),]

table(fwd$sample.id==rev$sample.id)
md<-cbind(fwd,rev$filename)
md$filename<-paste0(workdir,md$filename)
md$`rev$filename`<-paste0(workdir,md$`rev$filename`)
md<-md%>%
  arrange(sample.id)
colnames(md)<-c("sample-id",
                "animal",
                "forward-absolute-filepath",
                "reverse-absolute-filepath")
write.table(md,
            file = paste0(directory,"filenames-liu-local.tsv"),
            sep = "\t",row.names = F, quote = F)

####
directory<-"./data/bensch-data/metadata/"
workdir<-"/home/amir/bensch/data/raw/"
metadata<-read.table(paste0(directory,
                            "filenames-bensch.tsv"),
                     header = T,
                     sep = "\t")

fwd<-metadata[seq(from=1, to=nrow(metadata),by=2),]
rev<-metadata[seq(from=2, to=nrow(metadata),by=2),]

table(fwd$sample.id==rev$sample.id)
md<-cbind(fwd,rev$filename)
md$filename<-paste0(workdir,md$filename)
md$`rev$filename`<-paste0(workdir,md$`rev$filename`)
colnames(md)<-c("sample-id",
                "animal",
                "forward-absolute-filepath",
                "reverse-absolute-filepath")
write.table(md,
            file = paste0(directory,"filenames-bensch-raw-local.tsv"),
            sep = "\t",row.names = F, quote = F)
