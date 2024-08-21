#!/usr/bin/env bash
set -uex

FILE=raw_data/lambda_phage.fa
PARALLEL=0
IDRATE=0.1
FCOV=100

function sam2bam() {
    samtools sort -@20 "${1}".sam -o "${1}".bam
    samtools index "${1}".bam
    echo "${1}".bam
    python test_sam.py "${FILE}" "${1}".bam
}
mkdir -p tmp/
build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --i-file "${FILE}" \
    --read_len 125 \
    --mode template \
    --lc se \
    --i-parser memory \
    --i-fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq tmp/test_small_se_template_memory.fastq \
    --o-sam tmp/test_small_se_template_memory.sam \
    --o-hl_sam tmp/test_small_se_template.hl_memory.bam \
    --o-hl_sam-write_bam
sam2bam tmp/test_small_se_template_memory

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --i-file raw_data/ce11.mRNA.fa \
    --read_len 125 \
    --mode template \
    --lc pe \
    --i-parser memory \
    --i-fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_2 "${IDRATE}" \
    --o-fastq tmp/test_small_se_template_memory.fastq \
    --o-sam tmp/test_small_se_template_memory.sam
sam2bam tmp/test_small_se_template_memory
