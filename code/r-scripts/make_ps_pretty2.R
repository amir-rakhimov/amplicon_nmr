make_ps_pretty <- function(ps.object, agglom.rank) {
  ps.object[is.na(ps.object)] <- "NonNA"
  phylum <- which(colnames(ps.object) == "Phylum")
  agglom.rank.col <- which(colnames(ps.object) == agglom.rank)
  taxcols <- phylum:agglom.rank.col
  
  for (i in 1:nrow(ps.object)) {
    print(i)
    unclassified_indices <- which(ps.object[i, taxcols] %in% c("NonNA", "uncultured","uncultured_bacterium",
                                                               "unidentified_rumen","unidentified"))
    if (length(unclassified_indices) > 0) {
      previous_taxon <- ps.object[i, taxcols[min(unclassified_indices) - 1]]
      ps.object[i, taxcols[unclassified_indices]] <- paste0("Unclassified (", previous_taxon, ")")
    }
  }
  
  unclassified_cols <- ps.object[, agglom.rank.col] %in% c("Unclassified", "Uncultured")
  ps.object$Taxon <- ifelse(unclassified_cols, ps.object[[agglom.rank.col]], 
                            paste0(ps.object[[agglom.rank.col]], " (", ps.object[[agglom.rank.col - 1]], ")"))
  
  return(ps.object)
}
