library(Maaslin2)
?Maaslin2

# input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
# input_data
# 
# input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
# input_metadata
# 
# #get the pathway (functional) data - place holder
# download.file("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")

input_data = file.path("./data/maaslin2-tutorial","HMP2_taxonomy.tsv")
input_data

input_metadata = file.path("./data/maaslin2-tutorial","HMP2_metadata.tsv")
input_metadata

# HMP2_taxonomy.tsv: is a tab-demilited file with species as columns and samples 
# as rows. It is a subset of the taxonomy file so it just includes the species 
# abundances for all samples.
# 
# HMP2_metadata.tsv: is a tab-delimited file with samples as rows and metadata as
# columns. It is a subset of the metadata file so that it just includes some of 
# the fields.

# pathabundance_relab.tsv: is a tab-delimited file with samples as columns and 
# pathway level features as rows. It is a subset of the pathway file so it just 
# includes the species abundances for all samples.

# Saving inputs as data frames ####

df_input_data = read.table(file             = input_data,
                           header           = TRUE,
                           sep              = "\t", 
                           row.names        = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]

df_input_metadata = read.table(file             = input_metadata, 
                               header           = TRUE, 
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

df_input_path = read.csv("./data/maaslin2-tutorial/pathabundance_relab.tsv", 
                         sep              = "\t", 
                         stringsAsFactors = FALSE, 
                         row.names        = 1)
df_input_path[1:5, 1:5]

# 3.2 Running MaAsLin 2 ####
fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata, 
                    min_prevalence = 0,
                    # normalization  = "NONE", 
                    normalization = "TSS",
                    output         = "./data/maaslin2-tutorial/demo_output", 
                    fixed_effects  = c("diagnosis", "dysbiosis"),
                    reference      = c("diagnosis,nonIBD"),
                    transform = "AST",
                    analysis_method = "LM",
                    random_effects = NULL,
                    standardize = FALSE)

fit_data2 = Maaslin2(input_data     = df_input_data, 
                     input_metadata = df_input_metadata, 
                     min_prevalence = 0,
                     normalization  = "NONE",
                     output         = "./data/maaslin2-tutorial/demo_output2", 
                     fixed_effects  = c("diagnosis", "dysbiosis"),
                     reference      = c("diagnosis,nonIBD"))
