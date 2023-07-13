#!/bin/bash

#$ -V
#$ -cwd
#$ -l short
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -l d_rt=1:00:00
#$ -l s_rt=1:00:00
#$ -S /bin/bash

./prep_silva_data.py --infile song/SILVA_132_SSURef_tax_silva.fasta --taxafile song/silva-taxonomy.txt --outfasta silva-sequences.fasta --threads 10