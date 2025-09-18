#!/usr/bin/env bash
# shellcheck disable=SC2086
set -ue

SRC_DIR="$(dirname "$(readlink -f "${0}")")/../../"
cd "${SRC_DIR}" || exit 1
ORIG_TGZ_DIR="packing/artifacts/build_orig_tgz"
rm -fr "${ORIG_TGZ_DIR}"
mkdir -p "${ORIG_TGZ_DIR}"

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
        >"${ORIG_TGZ_DIR}"/git-files.txt

export BUILD_DIR="${ORIG_TGZ_DIR}/art-modern-${PACKAGE_VERSION}+dfsg"
mkdir -p "${BUILD_DIR}"

while read -r line; do
    if [ -f "${line}" ]; then
        install -D "${line}" "${BUILD_DIR}/${line}"
    fi
done <"${ORIG_TGZ_DIR}"/git-files.txt

env -C "${ORIG_TGZ_DIR}"/ \
    tar -czf "art-modern_${PACKAGE_VERSION}+dfsg.orig.tar.gz" \
    "art-modern-${PACKAGE_VERSION}+dfsg"
rm -fr "art-modern-${PACKAGE_VERSION}+dfsg"
