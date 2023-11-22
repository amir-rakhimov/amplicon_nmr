make_ps_pretty <- function(ps.object, agglom.rank) {
  ps.object[is.na(ps.object)] <- "NonNA"
  # find the column with your rank of interest
  # then extract columns between phylum and your rank
  phylum <- which(colnames(ps.object) == "Phylum") # identify the position of phylum column
  if (agglom.rank == "OTU") { # we don't trust species-level information
    agglom.rank.col <- which(colnames(ps.object) == "Genus")
  } else {
    agglom.rank.col <- which(colnames(ps.object) == agglom.rank)
  }
  taxcols <- phylum:agglom.rank.col # vector of taxonomic ranks
  # Replace empty taxa with Unclassified+previous taxon
  # we go through each taxonomic rank from phylum to agglom.rank
  # we substitute "uncultured" or "unclassified" and make it like 
  # Unclassified (Previous taxonomic rank)
  for (i in 1:nrow(ps.object)) {
    unclassified_indices <- which(ps.object[i, taxcols] %in% c("NonNA"))
    uncultured_indices <- which(ps.object[i, taxcols] %in% c("uncultured"))
    # we find the ranks that are "unclassified" or "uncultured"
    if (length(unclassified_indices) > 0 && length(uncultured_indices) == 0) {
      # if there is any unclassified, we find the name of previous taxonomic rank
      previous_taxon <- ps.object[i, taxcols[min(unclassified_indices) - 1]]
      if (length(previous_taxon)==0){
        # if all taxonomic ranks are unclassified except Kingdom, we set
        # all names to Unclassified (Bacteria)
        ps.object[i, taxcols[unclassified_indices]] <- paste0("Unclassified (", ps.object[i, phylum - 1], ")")
      }else{
        # otherwise, we use previous taxon name
        ps.object[i, taxcols[unclassified_indices]] <- paste0("Unclassified (", previous_taxon, ")")
      }
    }
    if (length(uncultured_indices) > 0 && length(unclassified_indices) == 0) {
      # if there is any uncultured, we find the name of previous taxonomic rank
      previous_taxon <- ps.object[i, taxcols[min(uncultured_indices) - 1]]
      if (length(previous_taxon)==0){
        # if all taxonomic ranks are Uncultured except Kingdom, we set
        # all names to Uncultured (Bacteria)
        ps.object[i, taxcols[uncultured_indices]] <- paste0("Uncultured (", ps.object[i, phylum - 1], ")")
      }else{
        # otherwise, we use previous taxon name
        ps.object[i, taxcols[uncultured_indices]] <- paste0("Uncultured (", previous_taxon, ")")
      }
    }
    
    ps.object[i,taxcols]<-gsub("\\_"," ", ps.object[i,taxcols])
    # gsub will remove underscores
  }
  # Create a Taxon column where items are either Unclassified (previous rank) 
  # or agglom.rank (previous taxonomic rank)
  # e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
  if(agglom.rank=="OTU"){
    ps.object$Taxon<-
      ifelse(grepl("Unclassified|Uncultured",
                   ps.object[["Genus"]]),  
             ps.object[["Genus"]], 
             paste0(ps.object[["Genus"]]," (",ps.object[["Family"]],")"))
  }else{
    ps.object$Taxon<-
      ifelse(grepl("Unclassified|Uncultured",
                   ps.object[[agglom.rank.col]]),  
             ps.object[[agglom.rank.col]], 
             paste0(ps.object[[agglom.rank.col]]," (",ps.object[[agglom.rank.col-1]],")"))
  }
  # grepl finds which agglom.rank column items are uncultured or unclassified
  # if true, keeps Unclassified (previous rank)
  # if false (so, it's not Unclassified, actually has taxonomic classification), 
  # converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)
  if(agglom.rank=="OTU"){
    ps.object$Species<-ps.object$Genus # species columns is same as genera
  }
  return(ps.object)
}
