#!/usr/bin/env bash
set -ue
curl https://ftp.cngb.org/pub/CNSA/data1/CNP0000126/CNS0020115/CNX0023695/CNR0028307/AF14-7LB_clean_500_100x.2.fastq.gz |
    gzip -cdf >data/e_coli.fa

curl -o data/cngb_aspera_download.key -s ftp://ftp.cngb.org/pub/Tool/Aspera/aspera_download.key
ascp \
    -i data/cngb_aspera_download.key -P 33001 -T -k 1 -l 100m \
    "aspera_download@183.239.175.39:/pub/CNSA/data1/CNP0000126/CNS0020115/CNX0023695/CNR0028307/AF14-7LB_clean_500_100x.1.fastq.gz" \
    data/e_coli_SCNR0028307_1.fastq.gz
ascp \
    -i data/cngb_aspera_download.key -P 33001 -T -k 1 -l 100m \
    "aspera_download@183.239.175.39:/pub/CNSA/data1/CNP0000126/CNS0020115/CNX0023695/CNR0028307/AF14-7LB_clean_500_100x.2.fastq.gz" \
    data/e_coli_SCNR0028307_2.fastq.gz

# PRJNA563564

seqtk seq -VQ 64 <(zcat data/e_coli_SCNR0028307_1.fastq.gz) | pigz -cf9 >data/e_coli_SCNR0028307_p33_1.fastq.gz
seqtk seq -VQ 64 <(zcat data/e_coli_SCNR0028307_2.fastq.gz) | pigz -cf9 >data/e_coli_SCNR0028307_p33_2.fastq.gz

mv data/e_coli_SCNR0028307_1.fastq.gz data/e_coli_SCNR0028307_1.fastq.gz.deleted
mv data/e_coli_SCNR0028307_2.fastq.gz data/e_coli_SCNR0028307_2.fastq.gz.deleted

bwa index data/e_coli.fa
bwa mem -t 20 \
    data/e_coli.fa \
    <(zcat data/e_coli_SCNR0028307_p33_1.fastq.gz) \
    <(zcat data/e_coli_SCNR0028307_p33_2.fastq.gz) |
    samtools sort -@ 20 --write-index -o data/e_coli_SCNR0028307.bam
perl src/pirs/baseCalling_Matrix_calculator.pl \
    -i data/e_coli_SCNR0028307.bam \
    -o data/e_coli_HiSeq2K_pirs_bcm \
    -r data/e_coli.fa \
    -m 0
samtools view -@20 -F 256 -F 2048 -F 4 data/e_coli.bam -o data/e_coli_primary_alignments.bam
perl src/pirs/indelstat_sam_bam.pl \
    data/e_coli_primary_alignments.bam \
    data/e_coli_pirs_indelstat
perl ../../deps/ART_profiler_illumina/art_profiler_illumina \
    data/e_coli_art_ data fastq.gz
