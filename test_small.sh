#!/usr/bin/env bash
set -ue

function sam2bam() {
    samtools sort "${1}"/1.sam -o "${1}"/1.bam
    samtools index "${1}"/1.bam
}

build/art \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --seq_file raw_data/small_test.fa \
    --out_file_prefix tmp/test_small_se \
    --read_len 125 \
    --fcov 10 \
    --parallel -1
sam2bam tmp/test_small_se

build/art \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --seq_file raw_data/small_test.fa \
    --out_file_prefix tmp/test_small_se_amp \
    --read_len 125 \
    --is_amplicon \
    --fcov 10 \
    --parallel -1
sam2bam tmp/test_small_se_amp

build/art \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file raw_data/small_test.fa \
    --out_file_prefix tmp/test_small_pe \
    --read_len 125 \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 10 \
    --is_pe \
    --fcov 10 \
    --parallel -1
sam2bam tmp/test_small_pe

build/art \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file raw_data/small_test.fa \
    --out_file_prefix tmp/test_small_pe_amp \
    --read_len 125 \
    --is_pe \
    --is_amplicon \
    --fcov 10 \
    --parallel -1
sam2bam tmp/test_small_pe_amp

build/art \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file raw_data/small_test.fa \
    --out_file_prefix tmp/test_small_mp \
    --read_len 125 \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 10 \
    --is_mp \
    --fcov 10 \
    --parallel -1
sam2bam tmp/test_small_mp

build/art \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
    --seq_file raw_data/small_test.fa \
    --out_file_prefix tmp/test_small_mp_amp \
    --read_len 125 \
    --is_mp \
    --is_amplicon \
    --fcov 10 \
    --parallel -1
sam2bam tmp/test_small_mp_amp
