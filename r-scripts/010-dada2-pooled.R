library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())

y.fastq.dir <- "./data/yasuda-data/fastq" 
y.final.dir <- file.path(y.fastq.dir, "final")
list.files(y.final.dir)
truncparams<-c(313,229)

# Getting ready ####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq 
# and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(y.final.dir, pattern="_clipped_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(y.final.dir, pattern="_clipped_L001_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles ####
plotQualityProfile(fnFs[1:2])+
  scale_x_continuous(breaks = seq(0,350, by=10))+
  theme(axis.text.x = element_text(angle=45,hjust=1))
plotQualityProfile(fnRs[1:2])+
  scale_x_continuous(breaks = seq(0,350, by=10))+
  theme(axis.text.x = element_text(angle=45,hjust=1))

# Filter and trim ####
## Place filtered files in filtered/ subdirectory ####
filtFs <- file.path(y.fastq.dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(y.fastq.dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# The maxEE parameter sets the maximum number of “expected errors” allowed in a 
# read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncparams,
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
View(out[1:nrow(out),2]/out[1:nrow(out),1]*100)

plotQualityProfile(filtFs[1:2])+
  scale_x_continuous(breaks = seq(0,350, by=10))+
  theme(axis.text.x = element_text(angle=45,hjust=1))
plotQualityProfile(filtRs[1:2])+
  scale_x_continuous(breaks = seq(0,350, by=10))+
  theme(axis.text.x = element_text(angle=45,hjust=1))

# Learn the error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## Visualise estimated error rates ####
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# Infer sample composition (Denoising) ####
# figuring out if each sequence is more likely to be of biological origin or 
# more likely to be spurious
# We are now ready to apply the core sample inference algorithm to the 
# filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Pooling data
dadaFs.pool<-dada(filtFs, err=errF, multithread=TRUE,
                  pool = TRUE)
dadaRs.pool<-dada(filtRs, err=errF, multithread=TRUE,
                  pool = TRUE)

dadaFs[[1]]
dadaFs.pool[[1]]

# dadaFs[[1]]: The DADA2 algorithm inferred 473 true sequence variants from the 
# 15224 unique sequences in the first sample. 
# dadaFs.pool[[1]]: 1076 sequence variants from 15224 input uniq sequences

dadaRs[[1]]
dadaRs.pool[[1]]

# Merge paired reads ####
mergers.pooled <- mergePairs(dadaFs.pool, filtFs, dadaRs.pool, filtRs, 
                             verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers.pooled[[1]])

# this object holds a lot of information that may be the first place you'd want 
# to look if you want to start poking under the hood
class(mergers.pooled) # list
length(mergers.pooled) # 23 elements in this list, one for each of our samples
names(mergers.pooled) # the names() function gives us the name of each element of the list 

class(mergers.pooled$D27) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe

names(mergers.pooled$D27) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" 
# "nindel"    "prefer"    "accept"


# Construct an amplicon sequence variant table (ASV) table ####
seqtab.pooled <- makeSequenceTable(mergers.pooled)
dim(seqtab.pooled)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.pooled)))

# Remove chimeras ####
seqtab.pooled.nochim <- removeBimeraDenovo(seqtab.pooled, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.pooled.nochim)

# though we only lost 3467 sequences, we don't know if they held a lot in 
# terms of abundance, this is one quick way to look at that

sum(seqtab.pooled.nochim)/sum(seqtab.pooled)
# 0.9688282 # in this case we barely lost any in terms of abundance


# Track reads through the pipeline ####
# As a final check of our progress, we’ll look at the number of reads that 
# made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track.pool <- cbind(out, sapply(dadaFs.pool, getN), sapply(dadaRs.pool, getN), 
               sapply(mergers.pooled, getN), rowSums(seqtab.pooled.nochim),
               round(rowSums(seqtab.nochim)/out[,1]*100, 1))

# If processing a single sample, remove the sapply calls: e.g. 
# replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track.pool) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim", "final_perc_reads_retained")
rownames(track.pool) <- sample.names
head(track.pool)

# Assign taxonomy ####
taxa.pooled <- 
  assignTaxonomy(seqtab.pooled.nochim, 
                 "./data/MiSeq_SOP/tax/silva_nr_v132_train_set.fa.gz", 
                 multithread=FALSE)

# Add species level assignment
taxa.pooled <- addSpecies(taxa.pooled, "./data/MiSeq_SOP/tax/silva_species_assignment_v132.fa.gz")

# Inspect taxonomic assignments
taxa.pooled.print <- taxa.pooled # Removing sequence rownames for display only
rownames(taxa.pooled.print) <- NULL
head(taxa.pooled.print)

# Phyloseq ####
# Sample names
samples.out <- rownames(seqtab.pooled.nochim)
# Class
animal.class<-c(rep("NMR",13),rep("Control",4),rep("NMR",6))
samdf.pool <- data.frame(Sample.id=samples.out, animal.class=animal.class)
rownames(samdf.pool) <- samples.out

## Construct the phyloseq object directly from dada2 output ####
ps.pooled <- phyloseq(otu_table(seqtab.pooled.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf.pool), 
               tax_table(taxa.pooled))

## use short names for ASVs ####
dna <- Biostrings::DNAStringSet(taxa_names(ps.pooled))
names(dna) <- taxa_names(ps.pooled)
ps.pooled <- merge_phyloseq(ps.pooled, dna)
taxa_names(ps.pooled) <- paste0("ASV", seq(ntaxa(ps.pooled)))
ps.pooled

# Visualise alpha diversity
plot_richness(ps.pooled, x="Sample.id", measures=c("Shannon", "Simpson"), 
              color="animal.class")

# Beta diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.pool.prop <- transform_sample_counts(ps.pooled, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.pool.prop, method="NMDS", distance="bray")

plot_ordination(ps.pool.prop, ord.nmds.bray, color="animal.class",
                title="Bray NMDS")

# Bar plot
top20 <- names(sort(taxa_sums(ps.pooled), decreasing=TRUE))[1:20]
ps.pool.top20 <- transform_sample_counts(ps.pooled, function(OTU) OTU/sum(OTU))
ps.pool.top20 <- prune_taxa(top20, ps.pool.top20)
plot_bar(ps.pool.top20, x="Sample.id", fill="Genus") + 
  facet_wrap(~animal.class, scales="free_x")
