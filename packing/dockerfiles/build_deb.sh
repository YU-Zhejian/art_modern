#!/usr/bin/env bash
set -ue

ORIG_TGZ="/mnt/art-modern_${PACKAGE_VERSION}+dfsg.orig.tar.gz"
DEBIAN_DIR="/mnt/debian"
# Extract this file to /mnt/build_deb/art-modern-${PACKAGE_VERSION}+dfsg
BUILD_DIR="/mnt/build_deb/art-modern-${PACKAGE_VERSION}+dfsg"
mkdir -p /mnt/build_deb
tar -xzf "${ORIG_TGZ}" -C /mnt/build_deb
cp -r "${DEBIAN_DIR}" "${BUILD_DIR}/debian"

# Build the deb
echo "DEB: Build the deb"
env -C "${BUILD_DIR}" debuild -us -uc \
    > >(sed --unbuffered 's/^/DEBUILD-INFO: /' | tee opt/build_deb/debuild-info.log) \
    2> >(sed --unbuffered 's/^/DEBUILD-ERR: /' | tee opt/build_deb/debuild-err.log >&2)
env -C "${BUILD_DIR}" debuild -- clean &>/dev/null
env -C "${BUILD_DIR}" debc 2>&1 | sed 's/^/DEBCHECK: /' | tee opt/build_deb/debcheck.log
env -C opt/build_deb/ \
    lintian art-modern_${PACKAGE_VERSION}+dfsg-1_amd64.changes \
    -EIi --pedantic --tag-display-limit 0 2>&1 |
    sed 's/^/LINTIAN: /' |
    tee opt/build_deb/lintian.log

