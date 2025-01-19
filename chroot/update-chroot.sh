#!/usr/bin/env bash
set -ue
SRC_DIR="$(dirname "$(readlink -f "${0}")")/../"
CHROOT_DIR="${1}"
cd "${SRC_DIR}" || exit 1
rm -rf "${CHROOT_DIR}"
mkdir -p "${CHROOT_DIR}"

git ls-files | while read -r line; do
    if [ -f "${line}" ]; then
        install -D "${line}" "${CHROOT_DIR}/${line}"
    fi
done
