#!/usr/bin/env bash
# mamba create -n bioconda -c conda-forge -c bioconda bioconda-utils
set -uxeo pipefail

conda run -n bioconda \
    --no-capture-output --live-stream \
    bioconda-utils lint \
    ./recepies \
    config.yml \
    --cache ./.cache
conda run -n bioconda \
    --no-capture-output --live-stream \
    bioconda-utils build \
    ./recepies \
    --docker \
    --mulled-test
