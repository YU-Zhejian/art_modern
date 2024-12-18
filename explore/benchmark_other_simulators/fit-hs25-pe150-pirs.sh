#!/usr/bin/env bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR584/DRR584679/DRR584679_1.fastq.gz \
    -O data/DRR584679_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR584/DRR584679/DRR584679_2.fastq.gz \
    -O data/DRR584679_2.fastq.gz
wget http://ftp.ensemblgenomes.org/pub/fungi/release-60/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz \
    -O data/yeast.fa.gz
gunzip data/yeast.fa.gz

bwa index data/yeast.fa
bwa mem -t 20 \
    data/yeast.fa \
    <(zcat data/DRR584679_1.fastq.gz) \
    <(zcat data/DRR584679_2.fastq.gz) | \
    samtools sort -@ 20 --write-index -o data/yeast.bam
perl src/pirs/baseCalling_Matrix_calculator.pl \
    -i data/yeast.bam \
    -o data/yeast_baseCalling_Matrix \
    -r data/yeast.fa \
    -m 0
samtools view -@20 -F 256 -F 2048 -F 4 data/yeast.bam -o data/yeast_primary_alignments.bam
perl src/pirs/indelstat_sam_bam.pl \
    data/yeast_primary_alignments.bam \
    data/yeast_indelstat_sam_bam
