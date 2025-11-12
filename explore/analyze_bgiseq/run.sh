#!/usr/bin/env bash
# SE BGISEQ-500
fasterq-dump SRR35105525 --threads 16 -O ./ --progress
fastqc SRR35105525.fastq

# PE BGISEQ-500
fasterq-dump -3 SRR7811627 --threads 16 -O ./ --progress
fastqc SRR7811627_1.fastq SRR7811627_2.fastq

for i in 1 2; do
    curl https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/047/SRR30664047/SRR30664047_"${i}".fastq.gz |
        gunzip -c | head -n 400000 >SRR30664047_"${i}".fastq
done
fastqc SRR30664047_1.fastq SRR30664047_2.fastq
