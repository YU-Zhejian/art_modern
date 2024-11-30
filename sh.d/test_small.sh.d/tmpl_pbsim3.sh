# shell=bash

parser=memory
coverage=pbsim3
for lc in se pe mp; do
    "${ART}" \
        --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
        --qual_file_2 data/Illumina_profiles/HiSeq2500L125R2.txt \
        --i-file data/raw_data/ce11.mRNA_head.pbsim3.transcript \
        --read_len 125 \
        --i-type pbsim3_template \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".sam
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}" "data/raw_data/ce11.mRNA_head.fa"
done

parser=stream
for lc in se pe mp; do
    "${ART}" \
        --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
        --qual_file_2 data/Illumina_profiles/HiSeq2500L125R2.txt \
        --i-file data/raw_data/ce11.mRNA_head.pbsim3.transcript \
        --read_len 125 \
        --i-type pbsim3_template \
        --i-batch_size 100 \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-hl_sam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam
    sam2bam "${OUT_DIR}"/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl "data/raw_data/ce11.mRNA_head.fa"
done
rm -fr "${OUT_DIR}"/test_small_??_template_"${parser}"_"${coverage}".hl.sam
