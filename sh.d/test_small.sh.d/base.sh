# shell=bash

export PARALLEL=0
export IDRATE=0.1

function sam2bam() {
    samtools sort -@20 --write-index "${1}".sam -o "${1}".bam
    python test_sam.py "${2}" "${1}".bam
    rm -f "${1}".sam "${1}".bam "${1}".bam.csi "${1}".bam.bai
}

mkdir -p tmp/
