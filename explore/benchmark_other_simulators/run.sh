#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

ART_MODERN_THREADS=36

function run() {
    /bin/time -a -o time.tsv -f '%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c\t%C' "${@}"
}

printf 'WALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\tCMD\n' >time.tsv
run bin/wgsim -1 150 -2 150 -N 3342880 -d 300 -s 20 data/ce11.fa data/wgsim_1.fq data/wgsim_2.fq >data/wgsim.vcf
run bin/art_original --seqSys HS25 --in data/ce11.fa -f 10 --len 150 --mflen 300 --sdev 20 --noALN --paired --out data/art
run opt/art_modern_build/art_modern \
    --mode wgs --lc pe \
    --i-file data/ce11.fa --i-fcov 10 --read_len 150 \
    --o-fastq data/art_modern_wgs_memory.fastq \
    --qual_file_1 ../../art/Illumina_profiles/HiSeq2500L150R1.txt \
    --qual_file_2 ../../art/Illumina_profiles/HiSeq2500L150R2.txt \
    --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel ${ART_MODERN_THREADS}

run bin/wgsim -1 150 -2 150 -N 185248 -d 300 -s 20 data/ce11_est.fa data/wgsim_1.fq data/wgsim_2.fq >data/wgsim.vcf
run bin/art_original --seqSys HS25 --in data/ce11_est.fa -f 10 --len 150 --mflen 300 --sdev 20 --noALN --paired --out data/art
run opt/art_modern_build/art_modern \
    --mode trans --lc pe \
    --i-file data/ce11_est.fa --i-fcov 10 --read_len 150 --i-parser memory \
    --o-fastq data/art_modern_wgs_memory.fastq \
    --qual_file_1 ../../art/Illumina_profiles/HiSeq2500L150R1.txt \
    --qual_file_2 ../../art/Illumina_profiles/HiSeq2500L150R2.txt \
    --pe_frag_dist_mean 300 --pe_frag_dist_std_dev 20 --parallel ${ART_MODERN_THREADS}
