#!/usr/bin/env bash
set -ue
mkdir -p ./monthly
aws s3 sync --no-sign-request s3://anaconda-package-data/conda/monthly/ monthly
