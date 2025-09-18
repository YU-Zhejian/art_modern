#!/usr/bin/env bash
set -ue
# Get mounted paths
ORIG_TGZ="/mnt/art-modern_${PACKAGE_VERSION}+dfsg.orig.tar.gz"
DEBIAN_DIR="/mnt/debian"

# Assemble the build directory
cp "${ORIG_TGZ}" /mnt/build_deb
env -C /mnt/build_deb tar -xzf "${ORIG_TGZ}"
BUILD_DIR="/mnt/build_deb/art-modern-${PACKAGE_VERSION}+dfsg"
cp -r "${DEBIAN_DIR}" "${BUILD_DIR}/debian"

# Set locales
export LANG=C
export LC_ALL=C
export LC_CTYPE=C
export LC_NUMERIC=C
export LC_TIME=C
export LC_COLLATE=C
export LC_MONETARY=C
export LC_MESSAGES=C
export LC_PAPER=C
export LC_NAME=C
export LC_ADDRESS=C
export LC_TELEPHONE=C
export LC_MEASUREMENT=C
export LC_IDENTIFICATION=C

# Build the deb
echo "DEB: Build the deb"
env -C "${BUILD_DIR}" debuild -us -uc \
    > >(sed --unbuffered 's/^/DEBUILD-INFO: /' | tee /mnt/build_deb/debuild-info.log) \
    2> >(sed --unbuffered 's/^/DEBUILD-ERR: /' | tee /mnt/build_deb/debuild-err.log >&2)
env -C "${BUILD_DIR}" debuild -- clean &>/dev/null
env -C "${BUILD_DIR}" debc 2>&1 | sed 's/^/DEBCHECK: /' | tee /mnt/build_deb/debcheck.log
env -C /mnt/build_deb/ \
    lintian art-modern_${PACKAGE_VERSION}+dfsg-1_amd64.changes \
    -EIi --pedantic --tag-display-limit 0 2>&1 |
    sed 's/^/LINTIAN: /' |
    tee /mnt/build_deb/lintian.log
