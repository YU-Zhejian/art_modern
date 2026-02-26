#!/usr/bin/env bash
set -ue
ARCH="$(uname -m)"
mkdir -p artifacts/RELEASE

./sh.d/prepare-master-dir-for-alpine.sh
ALPINE_MASTER_DIR="$(pwd)/artifacts/build_alpine_master"
singularity exec \
    --bind "${ALPINE_MASTER_DIR}":/mnt/art_modern-master \
    dockerfiles/alpine-latest.sif \
    sh /build_alpine_tgz.sh
mv artifacts/build_alpine_master/opt/build_rel_with_dbg_alpine-"${ARCH}".tar.gz \
    artifacts/RELEASE/build_rel_with_dbg_alpine-"${ARCH}".tar.gz
