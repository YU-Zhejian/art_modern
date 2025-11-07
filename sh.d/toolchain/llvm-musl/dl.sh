#!/usr/bin/env bash

# shellcheck disable=SC1091
set -ue
. vars.sh

rm -fr "${SRC_DIR}"
mkdir -p "${SRC_DIR}"
wget -4 https://musl.libc.org/releases/musl-1.2.5.tar.gz -O "${SRC_DIR}/musl-1.2.5.tar.gz"
env -C "${SRC_DIR}" tar -xf musl-1.2.5.tar.gz

wget -4 https://github.com/llvm/llvm-project/releases/download/llvmorg-18.1.8/llvm-project-18.1.8.src.tar.xz -O "${SRC_DIR}/llvm-project-18.1.8.src.tar.xz"
env -C "${SRC_DIR}" tar -xf llvm-project-18.1.8.src.tar.xz

wget -4 "https://www.kernel.org/pub/linux/kernel/v6.x/linux-6.14.9.tar.xz" -O "${SRC_DIR}/linux-6.14.9.tar.xz"
env -C "${SRC_DIR}" tar -xf linux-6.14.9.tar.xz

wget -4 https://github.com/ronchaine/libexecinfo/archive/d15ae4ca2e5e3b56d9abb2c7de03cb38078db474.zip -O "${SRC_DIR}/libexecinfo.zip"
env -C "${SRC_DIR}" unzip libexecinfo.zip

wget -4 https://ftpmirror.gnu.org/gnu/coreutils/coreutils-8.28.tar.xz -O "${SRC_DIR}/coreutils-8.28.tar.xz"
env -C "${SRC_DIR}" tar -xf coreutils-8.28.tar.xz

wget -4 "https://ftpmirror.gnu.org/gnu/bash/bash-5.3.tar.gz" -O "${SRC_DIR}/bash-5.3.tar.gz"
env -C "${SRC_DIR}" tar -xf bash-5.3.tar.gz
