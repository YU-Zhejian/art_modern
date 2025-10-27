#!/usr/bin/env bash
set -uxeo pipefail

mkdir -p ./test_output

conda run -n bioconda \
    --no-capture-output --live-stream \
    conda build purge-all
conda run -n bioconda \
    --no-capture-output --live-stream \
    conda build \
    --no-anaconda-upload \
    --strict-verify \
    --output-folder ./test_output \
    --package-format 2 \
    -c bioconda \
    -c conda-forge \
    ./recepies/art_modern \
    ./recepies/art_modern-openmpi
