#!/usr/bin/env bash
set -ue
export PACKAGE_VERSION=1.1.6 # FIXME
./sh.d/prepare-orig-tgz-for-deb.sh
for name in debian-12 debian-13 ubuntu-2404; do
    rm -fr artifacts/build_deb-"${name}"
    mkdir -p artifacts/build_deb-"${name}"
    singularity exec \
        --bind "artifacts/build_orig_tgz/art-modern_${PACKAGE_VERSION}+dfsg.orig.tar.gz":"/mnt/art-modern_${PACKAGE_VERSION}+dfsg.orig.tar.gz" \
        --bind ../debian:/mnt/debian \
        --bind artifacts/build_deb-"${name}":/mnt/build_deb \
        dockerfiles/${name}.sif \
        bash /build_deb.sh
    rm -fr artifacts/build_deb-"${name}"/"art-modern_${PACKAGE_VERSION}+dfsg"
done
