# shell=bash

build/art_modern \
    --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
    --sep_flag \
    --i-file raw_data/ce11.mRNA_head.fa \
    --read_len 125 \
    --mode template \
    --lc se \
    --i-parser memory \
    --i-fcov 1 \
    --parallel "${PARALLEL}" \
    --ins_rate_1 "${IDRATE}" \
    --del_rate_1 "${IDRATE}" \
    --o-fastq tmp/test_small_se_template_memory_sep.fastq \
    --o-sam tmp/test_small_se_template_memory_sep.sam \
    --o-hl_sam tmp/test_small_se_template_memory_sep.hl.bam \
    --o-hl_sam-write_bam
rm -fr tmp/test_small_se_template_memory_sep.*
