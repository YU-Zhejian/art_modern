#!/usr/bin/env bash
set -uxeo pipefail

conda run -n bioconda \
    --no-capture-output --live-stream \
    bioconda-utils lint \
    ./recepies \
    config.yml \
    --packages art_modern \
    --cache ./.cache \
    --logfile bioconda-utils-lint.log
conda run -n bioconda \
    --no-capture-output --live-stream \
    bioconda-utils build \
    ./recepies \
    --packages art_modern \
    --docker \
    --mulled-test \
    --logfile bioconda-utils-build.log
