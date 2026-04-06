#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -N 1
# Default mem 4GiB
#SBATCH --time=24:00:00
#SBATCH --job-name=raw_data-dl
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

module purge
module load petagene

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

set -ueo pipefail

mkdir -p raw_data
env -C raw_data pixi run -e moremodels axel \
    https://element-public-data.s3.us-west-2.amazonaws.com/20231031-2x300/human_wgs/DVT-0070/Samples/GAT-GEM-C001/GAT-GEM-C001_R1.fastq.gz
env -C raw_data pixi run -e moremodels axel \
    https://element-public-data.s3.us-west-2.amazonaws.com/20231031-2x300/human_wgs/DVT-0070/Samples/GAT-GEM-C001/GAT-GEM-C001_R2.fastq.gz

echo "DONE"
