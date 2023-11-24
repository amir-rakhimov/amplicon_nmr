library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())

path <- "./data/MiSeq_SOP/fastq" 
list.files(path)

# Getting ready ####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq 
# and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles ####
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim ####
## Place filtered files in filtered/ subdirectory ####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# The maxEE parameter sets the maximum number of “expected errors” allowed in a 
# read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
# On Windows set multithread=FALSE
# truncQ: (Optional). Default 2. Truncate reads at the first instance of a 
# quality score less than or equal to truncQ

# truncLen: (Optional). Default 0 (no truncation). Truncate reads after truncLen 
# bases.  Reads shorter than this are discarded

# maxN: (Optional). Default 0. After truncation, sequences with more than maxN 
# Ns will be discarded. Note that dada does not allow Ns.

# maxEE: (Optional). Default Inf (no EE filtering). After truncation, reads with 
# higher than maxEE "expected errors" will be discarded. Expected errors are 
# calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))

head(out)

# Learn the error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## Visualise estimated error rates ####
plotErrors(errF, nominalQ=TRUE)


# Sample inference ####
# We are now ready to apply the core sample inference algorithm to the 
# filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

# The DADA2 algorithm inferred 128 true sequence variants from the 1979 
# unique sequences in the first sample. 

# Merge paired reads ####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Construct an amplicon sequence variant table (ASV) table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline ####
# As a final check of our progress, we’ll look at the number of reads that 
# made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. 
# replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy ####
taxa <- assignTaxonomy(seqtab.nochim, "./data/MiSeq_SOP/tax/silva_nr_v132_train_set.fa.gz", 
                       multithread=FALSE)

# Add species level assignment
taxa <- addSpecies(taxa, "./data/MiSeq_SOP/tax/silva_species_assignment_v132.fa.gz")

# Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Evaluate accuracy ####
## Evaluating DADA2’s accuracy on the mock community: ####
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

# Phyloseq ####
# Sample names
samples.out <- rownames(seqtab.nochim)
# Subject names (F3 and Mock)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
# Gender
gender <- substr(subject,1,1)
# Subject id
subject <- substr(subject,2,999)
# day
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

## Construct the phyloseq object directly from dada2 output ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

## use short names for ASVs ####
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Visualise alpha diversity
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Beta diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

# Bar plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
