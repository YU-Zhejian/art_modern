#!/usr/bin/env bash
# shellcheck disable=SC2086
set -ue

SRC_DIR="$(dirname "$(readlink -f "${0}")")/../"
cd "${SRC_DIR}" || exit 1
rm -fr opt/build_deb
mkdir -p opt/build_deb

echo "DEB: Building for version ${PACKAGE_VERSION}"

# Generate orig tarball from git files
echo "DEB: Generate orig tarball from git files"
git ls-files |
    grep -v '^.idea/' |
    grep -v '^.github/' |
    grep -v '^.vscode/' |
    grep -v '^bioconda/' |
    grep -v '^chroot/' |
    grep -v '^debian/' |
    grep -v '^deps/ART_profiler_illumina/' |
    grep -v '^deps/concurrentqueue/' |
    grep -v '^deps/labw_slim_htslib/' |
    grep -v '^deps/pcg-cpp-0.98/' |
    grep -v '^deps/slim_fmt/' |
    grep -v '^deps/slim_abseil/' |
    grep -v '^deps/thread-pool/' |
    grep -v '^explore/' |
    grep -v '^env/' |
    grep -v '^sh.d/' \
        >opt/build_deb/git-files.txt
export BUILD_DIR="opt/build_deb/art-modern-${PACKAGE_VERSION}+dfsg"
mkdir -p "${BUILD_DIR}"
while read -r line; do
    if [ -f "${line}" ]; then
        install -D "${line}" "${BUILD_DIR}/${line}"
    fi
done <opt/build_deb/git-files.txt
env -C opt/build_deb/ \
    tar -czf "art-modern_${PACKAGE_VERSION}+dfsg.orig.tar.gz" \
    "art-modern-${PACKAGE_VERSION}+dfsg"

# Copy debian files
echo "DEB: Copy debian files"
cp -r debian "${BUILD_DIR}/debian"

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
