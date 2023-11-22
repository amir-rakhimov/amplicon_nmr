make_features_maaslin<-function(features.list, featurecolumn){
  features.list[[featurecolumn]]<-gsub(" \\(","\\.\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\)$","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("^\\[","X\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("^(\\d)", "X\\1", features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\] ","\\.\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\(\\[","\\(\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\[","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("-","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub(" ","\\.",features.list[[featurecolumn]])
  return(features.list)
}