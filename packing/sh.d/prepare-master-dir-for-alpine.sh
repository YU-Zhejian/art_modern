#!/usr/bin/env bash
set -ue
export ALPINE_MASTER_DIR="packing/artifacts/build_alpine_master"

SRC_DIR="$(dirname "$(readlink -f "${0}")")/../../"
cd "${SRC_DIR}" || exit 1
rm -fr "${ALPINE_MASTER_DIR}"
mkdir -p "${ALPINE_MASTER_DIR}"

echo "ALPINE: Generate master directory from git files"

git ls-files | while read -r line; do
    if [ -f "${line}" ]; then
        install -D "${line}" "${ALPINE_MASTER_DIR}/${line}"
    fi
done
