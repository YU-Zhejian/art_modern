#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=test-am-large-contigs
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

NJOBS=32
FCOV=20

# Set locale to C
export LC_ALL=C
export LANG=C


printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >time.tsv

# Only test whether the SAM output can be generated from such large contigs
"$(which time)" -a -o time.tsv -f '\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' \
    ../../opt/build_release_install/bin/art_modern \
    --i-file generated/large_contigs.fa \
    --i-type fasta \
    --i-parser htslib \
    --i-fcov "${FCOV}" \
    --parallel "${NJOBS}" \
    --o-hl_sam generated/large_contigs.fa.bam \
    --o-hl_sam-num_threads "${NJOBS}" \
    --o-hl_sam-write_bam

mkdir -p tmp
samtools view "generated/large_contigs.fa.bam" | \
    cut -f 1 | \
    sort -S 10K -T tmp | \
    uniq -d >"generated/large_contigs.fa.duplicated.txt"
cmp "generated/large_contigs.fa.duplicated.txt" /dev/null

