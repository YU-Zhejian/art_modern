#!/usr/bin/env bash
# shellcheck disable=SC2086
set -ue

export ORIG_TGZ_DIR="packing/artifacts/build_orig_tgz"
export ORIG_TGZ_PACKAGE_DIR="${ORIG_TGZ_DIR}/${PACKAGE_FULL_NAME}"

SRC_DIR="$(dirname "$(readlink -f "${0}")")/../../"
cd "${SRC_DIR}" || exit 1
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
    grep -v '^deps/slim_fmt/' |
    grep -v '^deps/slim_abseil/' |
    grep -v '^deps/thread-pool/' |
    grep -v '^explore/' |
    grep -v '^env/' |
    grep -v '^sh.d/' \
        >"${ORIG_TGZ_DIR}"/git-files.txt

mkdir -p "${ORIG_TGZ_PACKAGE_DIR}"

while read -r line; do
    if [ -f "${line}" ]; then
        install -D "${line}" "${ORIG_TGZ_PACKAGE_DIR}/${line}"
    fi
done <"${ORIG_TGZ_DIR}"/git-files.txt

env -C "${ORIG_TGZ_DIR}"/ \
    tar -czf "${PACKAGE_FULL_NAME}.orig.tar.gz" "${PACKAGE_FULL_NAME}"
rm -fr "${PACKAGE_FULL_NAME}"
