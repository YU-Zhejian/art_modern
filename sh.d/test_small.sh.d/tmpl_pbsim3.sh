# shell=bash

parser=memory
coverage=pbsim3
for lc in se pe mp; do
    build/art_modern \
        --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
        --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
        --i-file raw_data/ce11.mRNA_head.pbsim3.transcript \
        --read_len 125 \
        --i-type pbsim3_template \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-sam tmp/test_small_"${lc}"_template_"${parser}"_"${coverage}".sam
    sam2bam tmp/test_small_"${lc}"_template_"${parser}"_"${coverage}" raw_data/ce11.mRNA_head.fa
done

parser=stream
for lc in se pe mp; do
    build/art_modern \
        --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
        --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
        --i-file raw_data/ce11.mRNA_head.pbsim3.transcript \
        --read_len 125 \
        --i-type pbsim3_template \
        --i-batch_size 100 \
        --mode template \
        --lc "${lc}" \
        --i-parser "${parser}" \
        --parallel "${PARALLEL}" \
        --ins_rate_1 "${IDRATE}" \
        --del_rate_1 "${IDRATE}" \
        --o-hl_sam tmp/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam
    sam2bam tmp/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl raw_data/ce11.mRNA_head.fa
done
rm -fr tmp/test_small_??_template_"${parser}"_"${coverage}".hl.sam
