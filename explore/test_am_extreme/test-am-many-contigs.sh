#!/usr/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --mem=16G
#SBATCH --time=5:00:00
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
    mpirun -bootstrap=slurm -n "${SLURM_NTASKS}" \
        ../../opt/build_release_install-mpi/bin/art_modern-mpi \
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
