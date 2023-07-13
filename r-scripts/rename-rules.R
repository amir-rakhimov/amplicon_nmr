library(tidyverse)
workdir<-
f1<-"C:/Users/Artlab2/Documents/amir/phd_code/data/filenames-shanmu.tsv"
tab1<-read.table(f1, header = FALSE)
f2<-"C:/Users/Artlab2/Documents/amir/phd_code/data/shanmuganandam-data/metadata/metadata-clean.tsv"
tab2<-read.table(f2, header = TRUE)

tab1$Run<- gsub("_pass_[1,2].fastq.gz","",tab1$V1)

merged<-tab1%>%
  left_join(tab2, by = "Run")
merged$filename<-gsub("pass_", "R", merged$V1)
merged$read<-gsub("SRR[0-9]*_","",merged$filename)
merged$fullname<-paste(merged$sample.id,merged$read,sep = "_")

# metadata
f1<-"C:/Users/Artlab2/Documents/amir/phd_code/data/bensch-data/metadata/bensch-frozen.csv"
tab1<-read.table(f1, header = TRUE,sep = ",")
# old filenames
f2<-"C:/Users/Artlab2/Documents/amir/phd_code/data/bensch-data/metadata/fnames.txt"
tab2<-read.table(f2, header = FALSE)

tab2$Run<-gsub("_pass_[1,2].fastq.gz","",tab2$V1)


merged<-tab1%>%
  right_join(tab2,by=c("SRA"="Run"))
merged$filename<-gsub("pass_", "R", merged$V1)
merged$filename<-paste0("F",merged$sample,"_",merged$filename)

merged$read<-gsub("SRR[0-9]*_","",merged$filename)
merged$fullname<-paste0("F",merged$sample,"_",merged$read)

merged$sample.id<-paste0("F",merged$sample)
merged$animal<-"FukomysDamarensis"

rename_rules<-merged[,c("V1", "fullname")]
write.table(rename_rules,file = "C:/Users/Artlab2/Documents/amir/phd_code/data/bensch-data/metadata/bensch-rules.tsv", 
            sep="\t", 
            row.names = FALSE,col.names = FALSE, quote = F)

merged.to.save<-merged%>%
  select(sample.id,animal, fullname)
colnames(merged.to.save)<-c("sample-id", "animal", "filename")
write.table(merged.to.save,file = "./data/filenames-bensch.tsv", sep="\t", 
            row.names = FALSE,quote = F)
