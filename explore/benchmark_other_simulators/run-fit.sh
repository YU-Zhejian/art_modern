#!/usr/bin/env bash
set -ue

curl -o data/cngb_aspera_download.key -s ftp://ftp.cngb.org/pub/Tool/Aspera/aspera_download.key
for i in 1 2; do
    ascp \
        -i data/cngb_aspera_download.key -P 33001 -T -k 1 -l 100m \
        "aspera_download@183.239.175.39:/pub/CNSA/data1/CNP0000126/CNS0020115/CNX0023695/CNR0028307/AF14-7LB_clean_500_100x.${i}.fastq.gz" \
        data/e_coli_CNR0028307_"${i}".fastq.gz
    seqtk seq -VQ 64 <(zcat data/e_coli_CNR0028307_"${i}".fastq.gz) |
        pigz -cf9 \
            >data/e_coli_CNR0028307_p33_"${i}".fastq.gz
    mv data/e_coli_CNR0028307_"${i}".fastq.gz data/e_coli_CNR0028307_"${i}".fastq.gz.deleted

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR160/089/SRR16074289/SRR16074289_"${i}".fastq.gz \
        -O data/soybean_SRR16074289_"${i}".fastq.gz
done

mkdir -p data/e_coli_HiSeq2K data/soybean_HiSeq2500
mv data/e_coli_CNR0028307_p33_?.fastq.gz data/e_coli_HiSeq2K/
mv data/soybean_SRR16074289_?.fastq.gz data/soybean_HiSeq2500/

bash ../../deps/ART_profiler_illumina/art_profiler_illumina \
    data/e_coli_HiSeq2K_art_ data/e_coli_HiSeq2K fastq.gz

bash ../../deps/ART_profiler_illumina/art_profiler_illumina \
    data/soybean_HiSeq2500_art_ data/soybean_HiSeq2500 fastq.gz

bwa index data/e_coli.fa
bwa index data/soybean.fa

bwa mem -t 20 \
    data/e_coli.fa \
    <(zcat data/e_coli_HiSeq2K/e_coli_CNR0028307_p33_1.fastq.gz) \
    <(zcat data/e_coli_HiSeq2K/e_coli_CNR0028307_p33_2.fastq.gz) |
    samtools sort -@ 20 --write-index -o data/e_coli_CNR0028307.bam
perl src/pirs/baseCalling_Matrix_calculator.pl \
    -i data/e_coli_CNR0028307.bam \
    -o data/e_coli_HiSeq2K_pirs_bcm \
    -r data/e_coli.fa \
    -m 0
samtools view -@20 -F 256 -F 2048 -F 4 data/e_coli_CNR0028307.bam -o data/e_coli_CNR0028307_primary_alignments.bam
perl src/pirs/indelstat_sam_bam.pl \
    data/e_coli_CNR0028307_primary_alignments.bam \
    data/e_coli_HiSeq2K_pirs_indelstat

bwa mem -t 20 \
    data/soybean.fa \
    <(zcat data/soybean_HiSeq2500/soybean_SRR16074289_1.fastq.gz) \
    <(zcat data/soybean_HiSeq2500/soybean_SRR16074289_2.fastq.gz) |
    samtools sort -@ 20 --write-index -o data/soybean_SRR16074289.bam
perl src/pirs/baseCalling_Matrix_calculator.pl \
    -i data/soybean_SRR16074289.bam \
    -o data/soybean_HiSeq2500_pirs_bcm \
    -r data/soybean.fa \
    -m 0
samtools view -@20 -F 256 -F 2048 -F 4 data/soybean_SRR16074289.bam \
    -o data/soybean_SRR16074289_primary_alignments.bam
perl src/pirs/indelstat_sam_bam.pl \
    data/soybean_SRR16074289_primary_alignments.bam \
    data/soybean_HiSeq2500_pirs_indelstat

mkdir -p profiles
python get_profile.py data/e_coli_HiSeq2K/e_coli_CNR0028307_p33_1.fastq.gz profiles/pe100_real_1.profile 100 &
python get_profile.py data/e_coli_HiSeq2K/e_coli_CNR0028307_p33_2.fastq.gz profiles/pe100_real_2.profile 100 &

python get_profile.py data/soybean_HiSeq2500/soybean_SRR16074289_1.fastq.gz profiles/pe300_real_1.profile 300 &
python get_profile.py data/soybean_HiSeq2500/soybean_SRR16074289_2.fastq.gz profiles/pe300_real_2.profile 300 &
wait
