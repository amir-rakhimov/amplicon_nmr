# 16S rRNA gene sequencing analysis of the naked mole-rat gut microbiota
Original publication: Diversity and stability of the gut microbiome of naked mole-rat (Heterocephalus glaber), the longest-lived rodent  
Amir Rakhimov, Noriko Yasuda-Yoshihara, Masanori Arita, Kazuhiro Okumura, Yoshimi Kawamura, Kaori Oka, Hiroshi Mori, Yuichi Wakabayashi, Yoshifumi Baba, Hideo Baba, Kyoko Miura  
bioRxiv 2026.02.16.704739; doi: https://doi.org/10.64898/2026.02.16.704739


## Description  
This repository contains scripts and metadata related to the analysis of 
16S rRNA gene sequencing data from naked mole-rats and other rodents. 
Data from naked mole-rats and mice was generated at Kumamoto University and
Chiba Cancer Center Research Institute; data
from other hosts was taken from previous publications.

The R scripts were converted into R markdown files and knitr files, which can
be found in the markdown/ directory. The final report is _main.html

Scripts and metadata related to whole metagenome sequencing and MAG assembly can
be found in the other repository: https://github.com/amir-rakhimov/metagenome/

## Installation  
QIIME2 was run in a conda environment, requiring the following packages:

Statistical analysis in R requires the following R and Bioconductor packages:


# Usage
## QIIME2 pipeline  
For running code on a computer cluster with SLURM, run `.slurm` files from `code/slurm-scripts/`. For running code locally, run `.sh` scripts
from `code/bash-scripts/`  
1. Run `bash/001-get-fastq-files.sh` to download FASTQ files from other studies  
2. Run `qiime2-import-datasets` to import sequencing runs separately and trim primers with `q2-cutadapt`  
3. Run `qiime2-test-trunc_len` to test different DADA2 truncation parameters with `q2-dada2`  
4. Run `qiime2-train-classifier` to train a Naive Bayes classifier on SILVA database with `q2-feature-classifier`  
5. 

## Data analysis in R   
Run `.R` scripts in `code/r-scripts/`  


## Input files


## Example and output 
All the output can be found in the report at markdown/_main.html


## Publications
Original publication: Diversity and stability of the gut microbiome of naked mole-rat (Heterocephalus glaber), the longest-lived rodent  
Amir Rakhimov, Noriko Yasuda-Yoshihara, Masanori Arita, Kazuhiro Okumura, Yoshimi Kawamura, Kaori Oka, Hiroshi Mori, Yuichi Wakabayashi, Yoshifumi Baba, Hideo Baba, Kyoko Miura  
bioRxiv 2026.02.16.704739; doi: https://doi.org/10.64898/2026.02.16.704739









Amplicon_nmr project
