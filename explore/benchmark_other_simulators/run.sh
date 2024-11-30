#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
mkdir -p data_out

ART_MODERN_THREADS=36

function run() {
    /bin/time -a -o time.tsv -f '%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c\t%C' "${@}"
}

printf 'WALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\tCMD\n' >time.tsv
run bin/wgsim \
    -1 150 -2 150 -N 3342880 -d 300 -s 20 \
    data/ce11.fa \
    data_out/ce11_wgsim_1.fq data_out/ce11_wgsim_2.fq >data_out/wgsim.vcf
run bin/art_original \
    --in data/ce11.fa --out data_out/ce11_art_ \
    -f 10 --len 150 --mflen 300 --sdev 20 --noALN --paired --seqSys HS25
run opt/art_modern_build/art_modern \
    --mode wgs --lc pe \
    --i-file data/ce11.fa --i-fcov 10 --read_len 150 \
    --o-fastq data_out/ce11_art_modern_wgs_memory.fastq \
    --qual_file_1 ../../data/Illumina_profiles/HiSeq2500L150R1.txt \
    --qual_file_2 ../../data/Illumina_profiles/HiSeq2500L150R2.txt \
    --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel ${ART_MODERN_THREADS}

run bin/wgsim -1 150 -2 150 -N 10010885680 -d 300 -s 20 \
    data/hg38_long_mrna.fa \
    data_out/hg38_long_mrna_wgsim_1.fq data_out/hg38_long_mrna_wgsim_2.fq >data_out/hg38_long_mrna_wgsim.vcf
run bin/art_original \
    --in data/hg38_long_mrna.fa --out data_out/hg38_long_mrna_art_ \
    -f 10 --len 150 --mflen 300 --sdev 20 --noALN --paired --seqSys HS25
run opt/art_modern_build/art_modern \
    --mode trans --lc pe \
    --i-file data/hg38_long_mrna.fa --i-fcov 10 --read_len 150 --i-parser memory \
    --o-fastq data_out/hg38_long_mrna_art_modern_wgs_memory.fastq \
    --qual_file_1 ../../data/Illumina_profiles/HiSeq2500L150R1.txt \
    --qual_file_2 ../../data/Illumina_profiles/HiSeq2500L150R2.txt \
    --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel ${ART_MODERN_THREADS}
