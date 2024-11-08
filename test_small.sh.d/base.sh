# shell=bash

export PARALLEL=0
export IDRATE=0.1

function sam2bam() {
    samtools sort -@20 "${1}".sam -o "${1}".bam
    samtools index "${1}".bam
    echo "${1}".bam
    python test_sam.py "${2}" "${1}".bam
}

mkdir -p tmp/

# Test all output is working
build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --sep_flag \
    --i-file raw_data/ce11.mRNA_head.fa \
    --read_len 125 \
    --mode template \
    --lc se \
    --i-parser memory \
    --i-fcov 50 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq tmp/test_small_se_template_memory_sep.fastq \
    --o-sam tmp/test_small_se_template_memory_sep.sam \
    --o-hl_sam tmp/test_small_se_template_memory_sep.hl.bam \
    --o-hl_sam-write_bam
