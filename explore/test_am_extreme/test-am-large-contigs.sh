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

printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >time.tsv

# Only test whether the SAM output can be generated from such large contigs
"$(which time)" -a -o time.tsv -f '\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' \
    pixi run -e prev art_modern \
    --i-file generated/large_contigs.fa.gz \
    --i-type fasta \
    --i-parser htslib \
    --i-fcov "${FCOV}" \
    --parallel "${NJOBS}" \
    --o-sam generated/large_contigs.fa.gz.bam \
    --o-sam-num_threads "${NJOBS}" \
    --o-sam-write_bam
"$(which time)" -a -o time.tsv -f '\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' \
    pixi run -e prev art_modern \
    --i-file generated/large_contigs.fa \
    --i-type fasta \
    --i-parser htslib \
    --i-fcov "${FCOV}" \
    --parallel "${NJOBS}" \
    --o-sam generated/large_contigs.fa.bam \
    --o-sam-num_threads "${NJOBS}" \
    --o-sam-write_bam

for fn in \
    generated/large_contigs.fa.gz.bam \
    generated/large_contigs.fa.bam \
    ; do
    samtools sort -@ "${NJOBS}" --write-index -o "${fn%.bam}.sorted.bam" "${fn}"
    samtools depth -@ "${NJOBS}" -aa "${fn%.bam}.sorted.bam" \
        xz -9 -T"${NJOBS}" -vvvv > "${fn%.bam}.depth.tsv.xz"
done
