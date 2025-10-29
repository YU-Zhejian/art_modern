#!/usr/bin/env bash
set -ue
mkdir -p artifacts/RELEASE
ARCH="$(dpkg --print-architecture)"
if [ -z "${PACKAGE_VERSION:-}" ]; then
    export PACKAGE_VERSION="$(git describe --tags --abbrev=0)"
fi
export PACKAGE_FULL_NAME="art-modern_${PACKAGE_VERSION}+dfsg"

./sh.d/prepare-orig-tgz-for-deb.sh
for name in debian-13 ubuntu-2404; do
    rm -fr artifacts/build_deb-"${name}"
    mkdir -p artifacts/build_deb-"${name}"
    singularity exec \
        --bind "artifacts/build_orig_tgz/${PACKAGE_FULL_NAME}.orig.tar.gz":"/mnt/${PACKAGE_FULL_NAME}.orig.tar.gz" \
        --bind ../debian:/mnt/debian \
        --bind artifacts/build_deb-"${name}":/mnt/build_deb \
        dockerfiles/${name}.sif \
        bash /build_deb.sh
    rm -fr artifacts/build_deb-"${name}"/"${PACKAGE_FULL_NAME}"
    mv -v \
        artifacts/build_deb-"${name}"/"${PACKAGE_FULL_NAME}-1_${ARCH}.deb" \
        artifacts/RELEASE/"${PACKAGE_FULL_NAME}-1_${ARCH}-${name}.deb"
done
