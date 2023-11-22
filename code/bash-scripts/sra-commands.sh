../../sratoolkit.3.0.0-ubuntu64/bin/prefetch --option-file ../../song/sra-list.txt 


../../sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --outdir ../../song/data --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ./*sra

../../../sratoolkit.3.0.0-ubuntu64/bin/prefetch --option-file ../sra-list.txt --output-directory ./

../../../sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --outdir ../raw --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ./*sra
