# The function get_n_uniq_taxa_per_host groups the dataframe of abundances
# by host and taxonomic rank (e.g. Genus), retains unique rows (unique taxa in 
# each host), then groups by class and counts observations. It returns the
# number of unique taxa (e.g. genera) per host. Works for any number of hosts.
get_n_uniq_taxa_per_host<-function(tax.df,tax.rank){
  n.taxa.per.host<-tax.df%>%
    group_by_at(c("class",tax.rank))%>%
    distinct(get(tax.rank))%>%
    group_by(class)%>%
    tally
  return(n.taxa.per.host)
}