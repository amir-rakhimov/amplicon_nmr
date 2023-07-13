make_features_pretty<-function(features.list, featurecolumn){
  features.list[[featurecolumn]]<-gsub("\\.\\."," (",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\.$",")",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("X\\.","[",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\._","]_",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\(\\.","([",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\.","-",features.list[[featurecolumn]])
  return(features.list)
}