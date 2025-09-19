#!/usr/bin/env bash
set -ue
mkdir -p artifacts/RELEASE
ARCH=$(dpkg --print-architecture)
if [ -z "${PACKAGE_VERSION:-}" ]; then
    export PACKAGE_VERSION="$(git describe --tags --abbrev=0)"
fi
./sh.d/prepare-master-dir-for-alpine.sh
ALPINE_MASTER_DIR=artifacts/build_alpine_master
singularity exec \
    --bind "${ALPINE_MASTER_DIR}":/mnt/art_modern-master \
    dockerfiles/alpine-latest.sif \
    sh /build_alpine_tgz.sh
mv artifacts/build_alpine_master/opt/build_rel_with_dbg_alpine-x86_64.tar.gz \
    artifacts/RELEASE/build_rel_with_dbg_alpine-x86_64.tar.gz
