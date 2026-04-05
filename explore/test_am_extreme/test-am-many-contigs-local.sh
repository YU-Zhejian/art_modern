#!/usr/bin/bash

export PATH="${HOME}/.pixi/bin:${PATH}"
eval "$(pixi shell-hook)"

NJOBS=32
FCOV=20

if [ ! -f generated/many_contigs_36491262.fa ]; then
    pixi run -e testextreme \
        ./c/bin/generate_many_contigs | head -n 36491262 >generated/many_contigs_36491262.fa
    # 0x0011667ff
fi

. /opt/intel/oneapi/setvars.sh
set -ueo pipefail

collect=memory-consumption
rm -fr vtune-"${collect}"
vtune \
    -collect="${collect}" \
    -knob mem-object-size-min-thres=4096 \
    -source-search-dir="../../" \
    -result-dir=vtune-"${collect}" -- \
    /home/yuzj/Documents/art_modern/opt/build_release_install/bin/art_modern \
    --i-file generated/many_contigs_36491262.fa \
    --i-type fasta \
    --mode template \
    --lc se \
    --i-parser stream \
    --i-fcov "${FCOV}" \
    --i-batch_size 8192 \
    --parallel "${NJOBS}" \
    --reporting_interval-job_executor 10 \
    --reporting_interval-job_pool 50
vtune-gui vtune-"${collect}"
# 1048576
