library(tidyverse)
# library(phyloseq)
# library(Maaslin2)
# library(vegan)
# library(ALDEx2)
# NOT workspace dates but signif table dates
# maaslin.dates<-c("agegroup0_5"="20240613_16_01_22",
#                  "agegroup5_10"="20240613_16_02_23",
#                  "agegroup10_15"="20240613_16_03_39")
maaslin.dates<-c("agegroup0_10"="20240621_17_53_06")
maaslin.dates<-c("F"="20240621_17_54_53")
ref.level<-"agegroup0_10"
ref.level<-"F"

comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"

custom.levels<-c("agegroup0_10",
                 "agegroup10_16")
custom.levels<-c("F",
                 "M")


truncationlvl<-"234"
agglom.rank<-"OTU"
authorname<-"pooled"

rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"

host<-"NMR"
# host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")


host.labels<-
  c("NMR" = "*Heterocephalus glaber*")

# host.labels<-
#   c("B6mouse" = "B6 mouse",
#     "MSMmouse" = "MSM/Ms mouse",
#     "FVBNmouse" = "FVB/N mouse")

# load(file.path("./output/rdafiles",paste(
#   date_time,
#   authorname,read.end.type,"qiime2",
#   truncationlvl,agglom.rank,
#   "phyloseq-workspace.RData",sep = "-")))

# custom.levels<-c("agegroup0_5",
#                  "agegroup5_10",
#                  "agegroup10_15")


maaslin.signif.features<-list()
# aldex.signif.features<-list()
# ancombc.signif.features<-list()
# for (ref.level in custom.levels){
#   lvl.maaslin.signif.features<-
#     read.table(file.path("./rtables/pooled",
#                          paste("maaslin2",rare.status,
#                                filter.status,host,agglom.rank,
#                                comparison,truncationlvl,ref.level,
#                                "signif.tsv",sep="-")),
#                header = T,sep = "\t")
#   lvl.aldex.signif.features<-
#     read.table(file.path("./rtables/pooled",
#                          paste("aldex2",rare.status,
#                                filter.status,host,agglom.rank,
#                                comparison,truncationlvl,ref.level,
#                                "signif.tsv",sep="-")),
#                header = T,sep = "\t")
#   lvl.ancombc.signif.features<-
#     read.table(file.path("./rtables/pooled",
#                          paste("ancombc",rare.status,
#                                filter.status,host,agglom.rank,
#                                comparison,truncationlvl,ref.level,
#                                "signif.tsv",sep="-")),
#                header = T,sep = "\t")
#   lvl.ancombc.signif.features<-lvl.ancombc.signif.features%>%
#     dplyr::select(-X.Intercept.)%>%
#     pivot_longer(!taxon_id,names_to = "class",values_to = "coef")
#   
#   maaslin.signif.features[[ref.level]]<-lvl.maaslin.signif.features
#   aldex.signif.features[[ref.level]]<-lvl.aldex.signif.features
#   ancombc.signif.features[[ref.level]]<-lvl.ancombc.signif.features
# }

if(agglom.rank=="OTU"){
  for (ref.level in custom.levels){
    lvl.maaslin.signif.features<-
      read.table(file.path("./output/rtables/pooled",
                           paste(maaslin.dates[ref.level],"maaslin2",
                                 host,rare.status,
                                 filter.status,agglom.rank,
                                 comparison,truncationlvl,
                                 paste(sort(custom.levels),collapse = '-'),
                                 "ref",ref.level,"signif.tsv",sep="-")),
                 header = T,sep = "\t")
    # lvl.aldex.signif.features<-
    #   read.table(file.path("./output/rtables/pooled",
    #                        paste("aldex2",rare.status,
    #                              filter.status,host,agglom.rank,
    #                              comparison,truncationlvl,ref.level,
    #                              "signif.tsv",sep="-")),
    #              header = T,sep = "\t")
    # lvl.ancombc.signif.features<-
    #   read.table(file.path("./output/rtables/pooled",
    #                        paste("ancombc",rare.status,
    #                              filter.status,host,agglom.rank,
    #                              comparison,truncationlvl,ref.level,
    #                              "signif.tsv",sep="-")),
    #              header = T,sep = "\t")
    # lvl.ancombc.signif.features<-lvl.ancombc.signif.features%>%
    #   dplyr::select(-X.Intercept.)%>%
    #   pivot_longer(!taxon_id,names_to = "class",values_to = "coef")
    # 
    maaslin.signif.features[[ref.level]]<-lvl.maaslin.signif.features
    # aldex.signif.features[[ref.level]]<-lvl.aldex.signif.features
    # ancombc.signif.features[[ref.level]]<-lvl.ancombc.signif.features 
}
  }else{
 
}


rm(lvl.maaslin.signif.features)
# rm(lvl.aldex.signif.features)
# rm(lvl.ancombc.signif.features)

maaslin.signif.features<-bind_rows(maaslin.signif.features,.id = comparison)
# aldex.signif.features<-bind_rows(aldex.signif.features,.id = comparison)
# ancombc.signif.features<-bind_rows(ancombc.signif.features,.id = comparison)

maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by_at(c(comparison,"feature"))%>%
  filter(n()==length(custom.levels)-1)%>% 
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))

### For now, let's save maaslin output
# write.table(maaslin.signif.decreased,
#             file=file.path("./output/rtables/pooled",
#                            paste(
#                              paste(format(Sys.time(),format="%Y%m%d"),
#                                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "significant-decreased-all-lvls",host,rare.status,
#                              filter.status,agglom.rank,
#                              comparison,truncationlvl,
#                              paste(sort(custom.levels),collapse = '-'),
#                              ".tsv",sep="-")),
#             row.names = F,sep = "\t",col.names = T)
# 








# aldex.neg.effect<-aldex.signif.features%>%
#   filter(effect<0)

# ancombc.signif.decreased<-ancombc.signif.features%>%
#   as_tibble()%>%
#   group_by_at(c(comparison,"feature"))%>%
#   summarize(count=sum(coef<0))%>% # coef should be negative in all other groups
#   # so we count how many coefs are negative per taxon_id
#   filter(count==length(custom.levels)-1)%>%
#   left_join(ancombc.signif.features,by=c("strain"="strain", "taxon_id"="taxon_id"))

# Reduce(intersect,list(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="B6mouse"],
#                       aldex.neg.effect$OTU[aldex.neg.effect$strain=="B6mouse"],
#                       ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="B6mouse"]))

# ALDEX shares no otu with other two datasets simultaneously
# but separately, yes
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="B6mouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="B6mouse"])
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="B6mouse"],
          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="B6mouse"])
intersect(ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="B6mouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="B6mouse"])

intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="MSMmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="MSMmouse"])
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="MSMmouse"],
          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="MSMmouse"])
intersect(ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="MSMmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="MSMmouse"])

intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="FVBNmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="FVBNmouse"])
intersect(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain=="FVBNmouse"],
          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="FVBNmouse"])
intersect(ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain=="FVBNmouse"],
          aldex.neg.effect$OTU[aldex.neg.effect$strain=="FVBNmouse"])

signif.all.groups<-list()
for (ref.level in custom.levels){
  signif.lvl<-
    Reduce(intersect,list(maaslin.signif.decreased$feature[maaslin.signif.decreased$strain==ref.level],
                          ancombc.signif.decreased$taxon_id[ancombc.signif.decreased$strain==ref.level]))
  # filter by mean relative abundance: >=0.5%
  signif.lvl.filtered<-ps.q.agg%>%
    group_by(OTU,class)%>%
    filter(OTU %in% signif.lvl)%>%
    dplyr::select(OTU, MeanRelativeAbundance)%>%
    arrange(-MeanRelativeAbundance)%>%
    distinct(OTU,.keep_all = T)%>%
    ungroup()%>%
    pull(OTU)%>%
    unique()
  signif.all.groups[[ref.level]]<-signif.lvl.filtered
  write.table(signif.lvl.filtered,
              file=file.path("./rtables/pooled",
                             paste("significant-features",rare.status,
                                   filter.status,host,agglom.rank,
                                   comparison,truncationlvl,
                                   ref.level,"signif.tsv",sep="-")),
              row.names = F,sep = "\t",col.names = F)
}
rm(signif.lvl)

# save(signif.all.groups,
#      file=paste0("./rdafiles/",
#                paste("significant-features-all-groups",
#                      rare.status,filter.status,host,agglom.rank,
#                      comparison,truncationlvl,
#                      sep="-"),".RData"))


####
library(Polychrome)
library(ggtext)
ps.q.df.preprocessed.date_time<-"20240524_13_58_11" # for NMR
# "20240426_21_44_30" Rdata file for Genus all hosts
# "20240524_13_58_11" Rdata file for NMR
if(agglom.rank=="OTU"){
  ps.q.df.preprocessed<-read.table(
    file.path("./output/rtables",authorname,paste0(
      paste(
        ps.q.df.preprocessed.date_time,
        paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
        paste(names(host.labels),collapse = '-'),sep = "-"),
      ".tsv")),
    header = T)
}else{
  ps.q.df.preprocessed<-read.table(
    file.path("./output/rtables",authorname,paste0(
      paste(
        ps.q.df.preprocessed.date_time,
        "ps.q.df.rare-nonfiltered",agglom.rank,
        paste(custom.levels,collapse = '-'),sep = "-"),
      ".tsv")),
    header = T)}

if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class=="NMR",Abundance!=0)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))
  
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(agegroup=cut(age, breaks =c(0,10,16),
                        right = FALSE))
}else if(host=="mice"){
  # select mice and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(Abundance!=0)
  ps.q.df.preprocessed$agegroup<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
}

if (comparison=="age"){
  pretty.level.names<-names(table(ps.q.df.preprocessed$agegroup))
  names(pretty.level.names)<-ps.q.df.preprocessed%>%
    ungroup()%>%
    distinct(agegroup)%>%
    arrange(agegroup) %>%
    mutate(agegroup = paste0("agegroup", agegroup))%>%
    mutate(agegroup = gsub("\\(|\\)|\\[|\\]","",agegroup))%>%
    mutate(agegroup = gsub("\\,","_",agegroup))%>%
    pull()
  custom.levels<-names(pretty.level.names)
  gg.labs.name<-"Age group"
  gg.title.groups<-"age groups"
  
}else if (comparison=="sex"){
  pretty.level.names<-
    c("F" = "Females",
      "M" = "Males")
  custom.levels<-names(pretty.level.names)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  gg.labs.name<-"Host sex"
  gg.title.groups<-"groups"
}else if(comparison=="strain"){
  pretty.level.names<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.level.names),custom.md$class)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  gg.labs.name<-"Strain"
  gg.title.groups<-"strains"
}

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900"))
names(custom.fill)<-custom.levels
swatch(custom.fill)


ps.q.agg.date_time<-"20240620_12_38_18"
ps.q.agg<-readRDS(file.path("./output/rdafiles",paste(
  paste(ps.q.agg.date_time,
        "phyloseq-qiime",authorname,agglom.rank,read.end.type,truncationlvl,
        "table.rds",sep="-"))))


input_data_date_time<-"20240613_21_16_08"
load(file.path("./output/rdafiles",paste(
  input_data_date_time,
  "diffabund-input-data",host,rare.status,filter.status,agglom.rank,
  comparison,truncationlvl,
  paste(sort(custom.levels),collapse = '-'),"workspace.RData",sep="-")))


pretty.asv.names.df<-ps.q.agg%>%
  ungroup()%>%
  filter(get(agglom.rank)%in%maaslin.signif.features$feature,class==host)%>%
  distinct(get(agglom.rank),.keep_all = T)%>%
  select(OTU,Genus)%>%
  mutate(Taxon=paste0("ASV from ","<i>",Genus,"</i> (p = ",
                      round(maaslin.signif.features$qval,digits=3),")"))
pretty.asv.names<-c(pretty.asv.names.df$Taxon)
names(pretty.asv.names)<-pretty.asv.names.df$OTU

ps.q.agg%>%
  filter(get(agglom.rank)%in%maaslin.signif.features$feature,class==host)%>%
  left_join(custom.md[,c("Sample","agegroup")],by="Sample")%>%
  group_by_at(c(agglom.rank,"agegroup"))%>%
  ggplot(aes(x=factor(agegroup,level=c("agegroup0_10","agegroup10_16")),
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_boxplot(show.legend = FALSE)+
  facet_wrap(~OTU,scales = "free_y",
             ncol = 2,
             labeller = as_labeller(pretty.asv.names))+
  theme_bw()+
  labs(x="",
       y="Relative abundance (%)")+
 scale_color_manual(breaks =pretty.level.names,
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=c("agegroup0_10","agegroup10_16"))+ # rename boxplot labels (x axis)
  scale_fill_manual(values = custom.fill)+
  theme(axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x= ggtext::element_markdown(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "none")+
  ggtitle(paste0("Relative abundance of differentially abundant ASVs \nin different naked mole-rat groups"))
for(image.format in c("png","tiff")){
  ggsave(paste0("./images/taxaboxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      ref.level,"specific-bacteria",host,comparison,
                      sep = "-"),".",image.format),
         plot=last_plot(),
         width = 4000,height = 2000,
         units = "px",dpi=300,device = image.format)
}
