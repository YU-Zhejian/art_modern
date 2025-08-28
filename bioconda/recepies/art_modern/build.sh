#!/usr/bin/env bash
set -euo pipefail

mkdir -p /tmp/build_release
env -C /tmp/build_release cmake \
    -Wdev -Wdeprecated --warn-uninitialized \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCEU_CM_SHOULD_USE_NATIVE=OFF \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
    -DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    "${SRC_DIR}"
cmake --build /tmp/build_release --parallel ${CPU_COUNT}
cmake --install /tmp/build_release
