# shell=bash

for coverage in stranded strandless; do
    parser=memory
    for lc in se pe mp; do
        build/art_modern \
            --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
            --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
            --i-file raw_data/ce11.mRNA_head.fa \
            --read_len 125 \
            --mode template \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
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
            --i-file raw_data/ce11.mRNA_head.fa \
            --read_len 125 \
            --i-batch_size 100 \
            --mode template \
            --lc "${lc}" \
            --i-parser "${parser}" \
            --i-fcov raw_data/ce11.mRNA_head.cov_"${coverage}".tsv \
            --parallel "${PARALLEL}" \
            --ins_rate_1 "${IDRATE}" \
            --del_rate_1 "${IDRATE}" \
            --o-hl_sam tmp/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl.sam
        sam2bam tmp/test_small_"${lc}"_template_"${parser}"_"${coverage}".hl raw_data/ce11.mRNA_head.fa
    done
    rm -fr tmp/test_small_??_template_"${parser}"_"${coverage}".hl.sam
done
