---
output:
  html_document: default
  pdf_document: default
---
# PhD project: Cross-species gut microbiota comparison with a focus on naked mole-rat longevity

In this project, I analyse fecal microbiota from different rodents. Data from naked mole-rats (Heterocephalus glaber, NMR; n=19) and Specific-Pathogen-Free mice (Mus musculus, SPF mouse, n=4) were newly obtained at Kumamoto University. Data from MSM/Ms and FVB/N mice were obtained at Chiba Cancer Center. Reads from Damaraland mole-rats (Fukomys damarensis, DMR; n=20) [Bensch et al. 2022 PeerJ], European brown hares (Lepus europaeus, hare; n=9) and European rabbits (Oryctolagus cuniculus, rabbit; n=12) [Shanmuganandam et al. 2020 PeerJ], Lesser blind mole-rat (Nannospalax leucodon, spalax; n=15) [Sibai et al. 2020 OMICS J Integr Biol], and Siberian flying squirrel (Pteromys volans orii, PVO; n=19) [Liu et al. 2020 Sci Rep] were downloaded from previous studies. 

The repository is organised in four directories:

1. `code`: the folder with all scripts, such as bash and R scripts
- `r-scripts`:
- `bash-scripts`
2. `data`: raw data, metadata, and QIIME2 output that is used for downstream
processing in R.
- `fastq`: FASTQ files from amplicon sequencing. Each subfolder corresponds
to an experiment
  - `biagi-fastq`: from wild naked mole-rats (Data provided by Elena Biagi)
- `metadata`: metadata from each experiment
- `qiime`: QZA files from QIIME2, such as. reference tables, trees, and 
taxonomies. Each subfolder corresponds to an experiment
  - `biagi-qiime`: QIIME2 output from wild naked mole-rats
  - `pooled-qiime`: QIIME2 output from all experiments combined. We analysed all
  the raw data together
3. `images`: plots produced in R
- `barplots`: barplots of relative abundances
- `diversity`: alpha and beta diversity plots, PCA plots
  - `alpha`: alpha diversity plots
  - `nmds`: beta diversity analysis results visualised with nMDS
  - `pca`: PCA plots
  - `pcoa`: beta diversity analysis results visualised with PCoA
- `taxaboxplots`: boxplots of individual taxa
4. `output`: intermediate and final output from R
- `diffabund`: differential abundance analysis results
- `picrust`: PICRUSt2 results
- `rdafiles`: rdafiles from phyloseq, ALDEx2, MaAsLin2, ANCOM-BC
- `rtables`: R-tables from phyloseq, ALDEx2, MaAsLin2, ANCOM-BC
  - merged
  - pooled


