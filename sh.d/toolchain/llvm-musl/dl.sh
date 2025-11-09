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

wget -4 http://deb.debian.org/debian/pool/main/u/util-linux/util-linux_2.41.orig.tar.xz -O "${SRC_DIR}/util-linux_2.41.orig.tar.xz"
env -C "${SRC_DIR}" tar -xf util-linux_2.41.orig.tar.xz

wget -4 https://ftp.gnu.org/gnu/grep/grep-3.12.tar.xz -O "${SRC_DIR}/grep-3.12.tar.xz"
env -C "${SRC_DIR}" tar -xf grep-3.12.tar.xz

wget -4 https://ftp.gnu.org/gnu/sed/sed-4.9.tar.xz -O "${SRC_DIR}/sed-4.9.tar.xz"
env -C "${SRC_DIR}" tar -xf sed-4.9.tar.xz

# Note: 4.4 won't compile with musl, so we use 4.4.1
wget -4 https://ftp.gnu.org/gnu/make/make-4.4.1.tar.gz -O "${SRC_DIR}/make-4.4.1.tar.gz"
env -C "${SRC_DIR}" tar -xf make-4.4.1.tar.gz

wget -4 https://github.com/ninja-build/ninja/archive/refs/tags/v1.13.1.zip -O "${SRC_DIR}/ninja-1.13.1.zip"
env -C "${SRC_DIR}" unzip ninja-1.13.1.zip

wget https://www.zlib.net/zlib-1.3.1.tar.xz -O "${SRC_DIR}/zlib-1.3.1.tar.xz"
env -C "${SRC_DIR}" tar -xf zlib-1.3.1.tar.xz

wget https://distfiles.dereferenced.org/pkgconf/pkgconf-1.1.0.tar.xz -O "${SRC_DIR}/pkgconf-1.1.0.tar.xz"
env -C "${SRC_DIR}" tar -xf pkgconf-1.1.0.tar.xz

# NOTE: Align with system Python version to avoid compatibility issues
wget https://www.python.org/ftp/python/3.12.3/Python-3.12.3.tar.xz -O "${SRC_DIR}/Python-3.12.3.tar.xz"
env -C "${SRC_DIR}" tar -xf Python-3.12.3.tar.xz

wget https://github.com/tukaani-project/xz/releases/download/v5.8.1/xz-5.8.1.tar.xz -O "${SRC_DIR}/xz-5.8.1.tar.xz"
env -C "${SRC_DIR}" tar -xf xz-5.8.1.tar.xz

wget https://github.com/libffi/libffi/releases/download/v3.5.2/libffi-3.5.2.tar.gz -O "${SRC_DIR}/libffi-3.5.2.tar.gz"
env -C "${SRC_DIR}" tar -xf libffi-3.5.2.tar.gz

wget https://cmake.org/files/v3.17/cmake-3.17.5.tar.gz \
    -O "${SRC_DIR}/cmake-3.17.5.tar.gz"
env -C "${SRC_DIR}" tar -xf cmake-3.17.5.tar.gz

wget https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz -O "${SRC_DIR}/bzip2-1.0.8.tar.gz"
env -C "${SRC_DIR}" tar -xf bzip2-1.0.8.tar.gz
