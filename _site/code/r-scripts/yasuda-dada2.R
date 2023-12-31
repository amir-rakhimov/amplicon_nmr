library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())

y.fastq.dir <- "./data/yasuda-data/fastq" 
y.clipped.dir <- file.path(y.fastq.dir, "clipped")
list.files(y.clipped.dir)
truncparams<-c(313,229)
# Getting ready ####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq 
# and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(y.clipped.dir, pattern="_clipped_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(y.clipped.dir, pattern="_clipped_L001_R2_001.fastq.gz", full.names = TRUE))
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


dadaFs[[1]]
dadaRs[[1]]

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
               sapply(mergers, getN), rowSums(seqtab.nochim),
               round(rowSums(seqtab.nochim)/out[,1]*100, 1))

# If processing a single sample, remove the sapply calls: e.g. 
# replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim", "final_perc_reads_retained")
rownames(track) <- sample.names
head(track)

# Assign taxonomy ####
taxa <- 
  assignTaxonomy(seqtab.nochim, 
                 "./data/taxonomic reference/silva_nr99_v138.1_train_set.fa.gz", 
                 multithread=TRUE)

# Add species level assignment
taxa <- 
  addSpecies(taxa, 
             "./data/taxonomic reference/silva_species_assignment_v138.1.fa.gz")

# Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

trunclvl<-paste(truncparams,collapse  = "-")
save.image(paste0("./rdafiles/yasuda-r-dada2-",trunclvl,".RData"))

# Phyloseq ####
# Sample names
samples.out <- rownames(seqtab.nochim)
# Class
animal.class<-c(rep("NMR",13),rep("Control",4),rep("NMR",6))
samdf <- data.frame(Sample.id=samples.out, animal.class=animal.class)
rownames(samdf) <- samples.out

## Construct the phyloseq object directly from dada2 output ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

## use short names for ASVs ####
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Visualise alpha diversity
plot_richness(ps, x="Sample.id", measures=c("Shannon", "Simpson"), color="animal.class")

# Beta diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="animal.class", title="Bray NMDS")

# Bar plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample.id", fill="Genus") + 
  facet_wrap(~animal.class, scales="free_x")

nmr.ps<-subset_samples(ps,animal.class=="NMR")
nmr.top20<-names(sort(taxa_sums(nmr.ps), decreasing=TRUE))[1:20]
nmr.ps.top20 <- transform_sample_counts(nmr.ps, function(OTU) OTU/sum(OTU))
nmr.ps.top20 <- prune_taxa(nmr.top20, nmr.ps.top20)
plot_bar(nmr.ps.top20, x="Sample.id", fill="Genus") + 
  facet_wrap(~animal.class, scales="free_x")
