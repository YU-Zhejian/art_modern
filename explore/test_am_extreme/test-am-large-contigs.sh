#!/usr/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --job-name=test-am-large-contigs
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail -x

NJOBS=32
FCOV=20

# Set locale to C
export LC_ALL=C
export LANG=C
export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/current/lib64/libpmi2.so

# Only test whether the SAM output can be generated from such large contigs
pixi run -e buildoncluster \
    ../../opt/build_release_install-mpi/bin/art_modern-mpi \
    --i-file generated/large_contigs.fa \
    --i-type fasta \
    --i-parser htslib \
    --i-fcov "${FCOV}" \
    --parallel "${NJOBS}" \
    --o-hl_sam generated/large_contigs.fa.bam \
    --o-hl_sam-num_threads "${NJOBS}" \
    --o-hl_sam-write_bam

exit

mkdir -p tmp
samtools view "generated/large_contigs.fa.bam" | \
    cut -f 1 | \
    sort -S 10K -T tmp | \
    uniq -d >"generated/large_contigs.fa.duplicated.txt"
cmp "generated/large_contigs.fa.duplicated.txt" /dev/null

