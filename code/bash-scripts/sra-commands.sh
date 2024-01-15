#run `nano .bashrc`

# add this to the bottom of the file export PATH=sratoolkit.3.0.0-ubuntu64/bin:$PATH 

sra_list=../../song/sra-list.txt 
fastq_dir=../../song/data/fastq/raw


prefetch --option-file $sra_list --output-directory ./

fastq-dump --outdir $fastq_dir --gzip --skip-technical  --readids --read-filter pass \
--dumpbase --split-3 --clip ./*sra
