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
    https://element-public-data.s3.us-west-2.amazonaws.com/20240621-CloudbreakUltraQ/bases2fastq/MAXQ-0607__0fd9a6c2029c4f4fb6bf981f8a02c60f/Samples/GAT-ULTRAQ-L013/GAT-ULTRAQ-L013_R1.fastq.gz
env -C raw_data pixi run -e moremodels axel \
    https://element-public-data.s3.us-west-2.amazonaws.com/20240621-CloudbreakUltraQ/bases2fastq/MAXQ-0607__0fd9a6c2029c4f4fb6bf981f8a02c60f/Samples/GAT-ULTRAQ-L013/GAT-ULTRAQ-L013_R2.fastq.gz

echo "DONE"
