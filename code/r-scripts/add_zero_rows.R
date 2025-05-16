## Function to add zero rows if a species doesn't exist in a sample ####
# 1. Go through each sample
# 2. Find which organisms in the taxa.vector are missing in the sample
# 3. Create a tibble (new row)
# In order to have column names that are stored as string we make use 
# of bang bang operator !! which forces the evaluation of it succeeding name
# We also need to use walrus := instead of = which are equivalent and 
# prompts you to supply name (as is the case with our variable name) on it LHS (left hand side)
# tibble(Sample = "sample", 
#        Species =taxa,
#        Abundance = 0,
#        RelativeAbundance=0,
#        !!agglom.rank:="Treponema")
# 4. Add the row at the end of the dataframe
# 5. Fill in NA values in the empty columns based on non-empty rows in the Sample column 

add_zero_rows<-function(taxa.vector,tax.df,tax.rank){
  tax.df<-tax.df%>%ungroup()
  for (sample.name in unique(tax.df$Sample)) {
    missing_taxa <- setdiff(taxa.vector, tax.df[[tax.rank]][tax.df$Sample == sample.name])
    if (length(missing_taxa) > 0) {
      for (taxa in missing_taxa) {
        if(tax.rank=="Species"){
          new_row <- tibble(Sample = sample.name, 
                            Species= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0,
                            Genus=unique(pull(tax.df[which(tax.df$Species==taxa),"Genus"])))
        }else{
          new_row <- tibble(Sample = sample.name, 
                            !!tax.rank:= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0)
        }
        
        tax.df <- tax.df %>% add_row(.before = nrow(df), !!!new_row)
      }
    }
  }
  # To fill the NA values in the empty columns based on non-empty rows in the Sample column 
  # if(length(vars.to.fill!=0)){
  #   tax.df<- tax.df %>%
  #     group_by(Sample) %>%
  #     fill(all_of(vars.to.fill),.direction = "down")
  # }
  return(tax.df)
}

# The function add_zero_rows adds rows with Abundance = 0 to a dataframe tax.df
# based on a vector taxa.vector in the taxonomic rank tax.rank. When we create
# a barplot using the ggplot.species function, we want to show all bars, even 
# those with 0 value. So, we use the function add_zero_rows for adding empty bars.
# 1. The function loops through each sample in the dataframe: for (sample.name in unique(tax.df$Sample)).
# 2. The function creates a missing_taxa vector: setdiff(taxa.vector, tax.df[[tax.rank]][tax.df$Sample == sample.name]).
# The rationale is this: find which rows belong to a sample sample.name: tax.df$Sample == sample.name.
# The result is a vector of TRUE/FALSE values.
# Then, extract all rows with taxa at the tax.rank column as a vector: tax.df[[tax.rank]].
# Then, keep only those taxa that are found in the sample.name sample (if
# tax.df$Sample == sample.name is TRUE, then we keep tax.df[[tax.rank]] ).
# Finally, the function checks if the taxa in the taxa.vector are found in the filtered vector.
# If some taxon isn't found, the function will add a row with 0 Abundance. This
# taxon will be in the missing_taxa vector
# 3. Once the function creates the missing_taxa vector, it checks its length.
# If the length >0, the function will proceed to creation of zero rows: if (length(missing_taxa) > 0) 
# 4. While looping through each taxon in the missing_taxa vector, the function
# creates a new row as a tibble with 5 columns: Sample is the sample.name that 
# it's looping through in step 1. 
# The Abundance and RelativeAbundance columns will have a 0 value.
# If the tax.rank is Species, the next column is "Species", and the value is 
# the missing taxon. Otherwise, the function creates a column according to tax.rank.
# If the tax.rank is Species, the Genus column will have a value from the tax.df
# (tax.df is filtered by Species column and Genus is kep, then pooled): 
# unique(pull(tax.df[which(tax.df$Species==taxa),"Genus"])))
# 5. The new_row is added to the tax.df using the add_row function at the end 
# of the dataframe. !!! is a splicing operator that injects a list of 
# arguments: tax.df %>% add_row(.before = nrow(df), !!!new_row)
# https://rlang.r-lib.org/reference/splice-operator.html
# https://rlang.r-lib.org/reference/topic-inject.html
# !!! is not the same as !!
# 6. The columns that have NA value after injection will be filled by the 
# corresponding value from the same sample. So, if "age" column is empty in the
# 2D10 sample, dplyr will take the value from the other row that is not empty.

