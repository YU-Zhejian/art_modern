#!/usr/bin/env bash
set -uex

CONTIG=NC_001416.1
FILE=raw_data/lambda_phage.fa
PARALLEL=0
IDRATE=0.1
FCOV=100

function sam2bam() {
    samtools sort -@20 "${1}"/"${CONTIG}".sam -o "${1}"/"${CONTIG}".bam
    samtools index "${1}"/"${CONTIG}".bam
    echo "${1}"/"${CONTIG}".bam
    python test_sam.py "${FILE}" "${1}"/"${CONTIG}".bam
}

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --seq_file "${FILE}" \
    --out_file_prefix tmp/test_small_se_amp \
    --read_len 125 \
    --mode template \
    --fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --parallel_on_read \
    --ins_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --del_rate_2 "${IDRATE}"
sam2bam tmp/test_small_se_amp

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --seq_file "${FILE}" \
    --out_file_prefix tmp/test_small_se \
    --read_len 125 \
    --fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --parallel_on_read \
    --ins_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --del_rate_2 "${IDRATE}"
sam2bam tmp/test_small_se
# Official
time art_illumina \
    -ss HS25 \
    -sam \
    -i raw_data/lambda_phage.fa \
    -l 125 \
    -f "${FCOV}" \
    -o tmp/test_small_se_official \
    -ir "${IDRATE}" -ir2 "${IDRATE}" -dr "${IDRATE}" -dr2 "${IDRATE}"

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file "${FILE}" \
    --out_file_prefix tmp/test_small_pe \
    --read_len 125 \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 10 \
    --lc pe \
    --fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --parallel_on_read \
    --ins_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --del_rate_2 "${IDRATE}"
sam2bam tmp/test_small_pe

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file "${FILE}" \
    --out_file_prefix tmp/test_small_pe_amp \
    --read_len 125 \
    --lc pe \
    --mode template \
    --fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --parallel_on_read \
    --ins_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --del_rate_2 "${IDRATE}"
sam2bam tmp/test_small_pe_amp

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file "${FILE}" \
    --out_file_prefix tmp/test_small_mp \
    --read_len 125 \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 10 \
    --lc mp \
    --fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --parallel_on_read \
    --ins_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --del_rate_2 "${IDRATE}"
sam2bam tmp/test_small_mp

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file "${FILE}" \
    --out_file_prefix tmp/test_small_mp_amp \
    --read_len 125 \
    --lc mp \
    --mode template \
    --fcov "${FCOV}" \
    --parallel "${PARALLEL}" \
    --parallel_on_read \
    --ins_rate_1 "${IDRATE}" \
    --ins_rate_2 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --del_rate_2 "${IDRATE}"
sam2bam tmp/test_small_mp_amp
