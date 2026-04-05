#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=test-am-many-contigs
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

# FIXME: OOM Killed.

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

NJOBS=32
FCOV=20

pixi run -e testextreme \
    ./c/bin/generate_many_contigs |
    pixi run -e prev art_modern \
        --i-file /dev/stdin \
        --i-type fasta \
        --mode template \
        --lc se \
        --i-parser stream \
        --i-fcov "${FCOV}" \
        --i-batch_size 1048576 \
        --parallel "${NJOBS}" \
        --o-hl_sam generated/many_contigs.bam \
        --o-hl_sam-num_threads "${NJOBS}" \
        --o-hl_sam-compress_level 9 \
        --o-hl_sam-write_bam \
        --reporting_interval-job_executor 10 \
        --reporting_interval-job_pool 50
samtools fasta generated/many_contigs.bam >generated/many_contigs.fa
# This generates 10G reads.
