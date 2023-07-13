make_ps_pretty<-function(ps.object,agglom.rank){
  ps.object[is.na(ps.object)]<-"NonNA"
  # find the column with your rank of interest
  # then extract columns between phylum and your rank
  phylum<-which(colnames(ps.object) =="Phylum") # identify the position of phylum column
  agglom.rank.col<-which(colnames(ps.object) ==agglom.rank) # identify the 
  # position of agglom.rank column
  taxcols<-phylum:agglom.rank.col # vector of taxonomic ranks
  
  # Replace empty taxa with Unclassified+previous taxon
  # we go through each taxonomic rank from phylum to agglom.rank
  # we substitute "uncultured" or "unclassified" and make it like 
  # Unclassified (Previous taxonomic rank)

  for (i in 1:nrow(ps.object)){
    for (j in taxcols){
      if(ps.object[i,j]=="NonNA"){
        ps.object[i,j]<-paste0("Unclassified"," (",ps.object[i,j-1],")") 
      }
      if (grepl("Unclassified",ps.object[i, j-1])){
        ps.object[i,j]<-ps.object[i,j-1]
      }
      if (ps.object[i,j]=="uncultured"){
        ps.object[i,j]<-paste0("Uncultured"," (",ps.object[i,j-1],")")
      }
      if (grepl("Uncultured",ps.object[i, j-1])){
        ps.object[i,j]<-ps.object[i,j-1]
      }
    }
  }
  # Create a Taxon column where items are either Unclassified (previous rank) 
  # or agglom.rank (previous taxonomic rank)
  # e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
  ps.object$Taxon<-
    ifelse(grepl("Unclassified|Uncultured",
                 ps.object[[agglom.rank.col]]),  
           ps.object[[agglom.rank.col]], 
           paste0(ps.object[[agglom.rank.col]]," (",ps.object[[agglom.rank.col-1]],")")) 
  # grepl finds which agglom.rank column items are uncultured or unclassified
  # if true, keeps Unclassified (previous rank)
  # if false (so, it's not Unclassified, actually has taxonomic classification), 
  # converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)
  return(ps.object)
  
}
