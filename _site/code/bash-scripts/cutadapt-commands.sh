# How to see a list of files that have "_L001_R1_001.fastq.gz" at the end
# sed will substitute the part after s/
for i in *_L001_R1_001.fastq.gz
do
  SAMPLE=$(echo ${i} | sed "s/_R1_001\.fastq\.gz//")
  echo ${SAMPLE}_L001_R1_001.fastq.gz ${SAMPLE}_L001_R2_001.fastq.gz
done

# How to see sample names without showing the "_L001_R1_001.fastq.gz" part
# sed will substitute the part after s/
for i in *_L001_R1_001.fastq.gz
do
  echo ${i} | sed "s/_L001_R1_001\.fastq\.gz//";
done

# Cutadapt: Paired end reads
## regular trimming
# two 5' adapters
for i in raw/*_L001_R1_001.fastq.gz
do
  sample=$(echo ${i} | sed "s/_L001_R1_001\.fastq\.gz//" | sed "s/raw\///")
  cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o clipped5/${sample}_clipped5_L001_R1_001.fastq.gz -p clipped5/${sample}_clipped5_L001_R2_001.fastq.gz raw/${sample}_L001_R1_001.fastq.gz raw/${sample}_L001_R2_001.fastq.gz;
done

# Trim 3' adapters
for i in clipped5/*_clipped5_L001_R1_001.fastq.gz
do
  sample=$(echo ${i} | sed "s/_clipped5_L001_R1_001\.fastq\.gz//" | sed "s/clipped5\///")
  cutadapt -a GGATTAGATACCCBDGTAGTC$ -A CTGCWGCCNCCCGTAGG$ -o final/${sample}_final_L001_R1_001.fastq.gz -p final/${sample}_final_L001_R2_001.fastq.gz clipped5/${sample}_clipped5_L001_R1_001.fastq.gz clipped5/${sample}_clipped5_L001_R2_001.fastq.gz;
done

############################
## anchored trimming
# trim 5'
for i in raw/*_L001_R1_001.fastq.gz
do
  sample=$(echo ${i} | sed "s/_L001_R1_001\.fastq\.gz//" | sed "s/raw\///")
  cutadapt -g ^CCTACGGGNGGCWGCAG -G ^GACTACHVGGGTATCTAATCC -o anchored/${sample}_anchored_L001_R1_001.fastq.gz -p anchored/${sample}_anchored_L001_R2_001.fastq.gz raw/${sample}_L001_R1_001.fastq.gz raw/${sample}_L001_R2_001.fastq.gz;
done

# Trim 3' adapters
for i in clipped5/*_clipped5_L001_R1_001.fastq.gz
do
  sample=$(echo ${i} | sed "s/_clipped5_L001_R1_001\.fastq\.gz//" | sed "s/clipped5\///")
  cutadapt -a GGATTAGATACCCBDGTAGTC$ -A CTGCWGCCNCCCGTAGG$ -o final/${sample}_final_L001_R1_001.fastq.gz -p final/${sample}_final_L001_R2_001.fastq.gz clipped5/${sample}_clipped5_L001_R1_001.fastq.gz clipped5/${sample}_clipped5_L001_R2_001.fastq.gz;
done
######################################
# Single end reads
# regular trimming: start from the 5' end. Output in the folder clippedfront/
for i in raw/*.fastq.gz
do
  sample=$(echo ${i} | sed "s/\.fastq\.gz//" | sed "s/raw\///");
  cutadapt --report=minimal -g CCTACGGGNGGCWGCAG -o clippedfront/${sample}_clipped.fastq.gz raw/${sample}.fastq.gz;
done

# 3' end trimming: use files in clippedfront/ and output into final/
# use reverse complement to trim the 3' end
for i in clippedfront/*.fastq.gz
do
  sample=$(echo ${i} | sed "s/\.fastq\.gz//" | sed "s/clippedfront\///");
  cutadapt --report=minimal -a GGATTAGATACCCBDGTAGTC -o final/${sample}_final.fastq.gz clippedfront/${sample}.fastq.gz;
done
######################
# This is just to see trim file
cutadapt --report=minimal -g CCTACGGGNGGCWGCAG -o T01_clipped.fastq.gz T01.fastq.gz
cutadapt --report=minimal -a GGATTAGATACCCBDGTAGTC -o T01_final.fastq.gz T01_clipped.fastq.gz

# Primers to search for when using less
Amplicon PCR Forward Primer 
5′-CCTACGGGNGGCWGCAG-3′ (CCTACGGG[ACGT]GGC[AT]GCAG)
Amplicon PCR Reverse Primer
5′-GACTACHVGGGTATCTAATCC-3′ (GACTAC[ACT][ACG]GGGTATCTAATCC)

forward rev comp
CTGCWGCCNCCCGTAGG (CTGC[AT]GCC[ACGT]CCCGTAGG)
reverse rev comp
GGATTAGATACCCBDGTAGTC (GGATTAGATACCC[CGT][AGT]GTAGTC)

-g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC
-a GGATTAGATACCCBDGTAGTC$ -A CTGCWGCCNCCCGTAGG$


