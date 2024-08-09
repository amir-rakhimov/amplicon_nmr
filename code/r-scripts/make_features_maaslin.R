make_features_maaslin<-function(features.list, featurecolumn){
  features.list[[featurecolumn]]<-gsub("\\(","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\)","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\:","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\)$","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("^\\[","X\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("^(\\d)", "X\\1", features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\]","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\(\\[","\\(\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\[","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\]$","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("-","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub("\\+","\\.",features.list[[featurecolumn]])
  features.list[[featurecolumn]]<-gsub(" ","\\.",features.list[[featurecolumn]])
  return(features.list)
}
#original DIHYDRODIPICSYN-RXN: (expasy) 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7]
# 1. DIHYDRODIPICSYN-RXN: .expasy) 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7]
# 2. DIHYDRODIPICSYN-RXN: .expasy. 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7]
# 3. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7]
# 4. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7] 
#nothing here but end of the line ) will become .
# 5. DIHYDRODIPICSYN-RXN. .expasy .4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7] 
# nothing here but start of the line [ will become X
# 6. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7] 
# I forgot
# 7. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7] 
# nothing here but ] with space becomes ..
# 8. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase [4.3.3.7] 
# nothing here but ([ becomes ..
# 9. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase .4.3.3.7] 
# 10. DIHYDRODIPICSYN-RXN. .expasy. 4-hydroxy-tetrahydrodipicolinate synthase .4.3.3.7.
# 11. DIHYDRODIPICSYN.RXN. .expasy. 4.hydroxy.tetrahydrodipicolinate synthase .4.3.3.7. 
# 12. DIHYDRODIPICSYN.RXN. .expasy. 4.hydroxy.tetrahydrodipicolinate synthase .4.3.3.7. 
# nothing here but + becomes .
# 13. DIHYDRODIPICSYN.RXN...expasy..4.hydroxy.tetrahydrodipicolinate.synthase..4.3.3.7. 
