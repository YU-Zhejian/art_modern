#!/usr/bin/env bash
set -ue
ARCH="$(uname -m)"
ARCH=$(dpkg --print-architecture)
if [ -z "${PACKAGE_VERSION:-}" ]; then
    export PACKAGE_VERSION="$(git describe --tags --abbrev=0)"
fi
./sh.d/prepare-master-dir-for-alpine.sh
ALPINE_MASTER_DIR=artifacts/build_alpine_master
singularity exec \
    --bind "${ALPINE_MASTER_DIR}":/mnt/art_modern-master \
    dockerfiles/alpine-latest.sif \
    sh -c 'OPT_DIR=opt/; cd /mnt/art_modern-master && mkdir -p ${OPT_DIR}/build_rel_with_dbg_alpine && \
                                        	env -C ${OPT_DIR}/build_rel_with_dbg_alpine cmake \
                                        		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
                                        		-DCEU_CM_SHOULD_ENABLE_TEST=OFF \
                                        		-DCEU_CM_SHOULD_USE_NATIVE=OFF \
                                        		-DCMAKE_VERBOSE_MAKEFILE=ON \
                                        		-DBUILD_SHARED_LIBS=OFF \
                                        		-DUSE_MALLOC=NOP \
                                        		-DCMAKE_INSTALL_LIBDIR=bin \
                                        		-DCMAKE_INSTALL_INCLUDEDIR=bin \
                                        		-DCMAKE_INSTALL_PREFIX=${OPT_DIR}/build_rel_with_dbg_alpine_install /mnt/art_modern-master; cmake --build ${OPT_DIR}/build_rel_with_dbg_alpine --target art_tsam2gsam -- -j $(nproc); cmake --install .'

