#!/usr/bin/env bash
set -uxeo pipefail

mkdir -p ./test_output

conda run -n bioconda \
    --no-capture-output --live-stream \
    conda build purge-all
conda run -n bioconda \
    --no-capture-output --live-stream \
    conda build ./recepies/art_modern \
    --no-anaconda-upload \
    --strict-verify \
    --output-folder ./test_output \
    --package-format 2 \
    -c bioconda \
    -c conda-forge
