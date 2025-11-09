# shellcheck shell=bash
# shellcheck disable=SC2034

PROJ_DIR="$(pwd)"
SRC_DIR="${PROJ_DIR}/src"
MUSL_SRC_DIR="${SRC_DIR}/musl-1.2.5"
LLVM_SRC_DIR="${SRC_DIR}/llvm-project-18.1.8.src"
LLVM_CMAKE_SRC_DIR="${LLVM_SRC_DIR}/llvm"
LLVM_COMPILER_RT_SRC_DIR="${LLVM_SRC_DIR}/compiler-rt"
LINUX_SRC_DIR="${SRC_DIR}/linux-6.14.9"
BUSYBOX_SRC_DIR="${SRC_DIR}/busybox-upstream-1.37.0"
COREUTILS_SRC_DIR="${SRC_DIR}/coreutils-8.28"
BASH_SRC_DIR="${SRC_DIR}/bash-5.3"
LIBEXECINFO_SRC_DIR="${SRC_DIR}/libexecinfo-d15ae4ca2e5e3b56d9abb2c7de03cb38078db474"

SYSROOT="${PROJ_DIR}/sysroot"
KERNEL_PARAMS="nokaslr console=tty0 console=ttyS0,115200"

HOSTCC="$(which clang)"
HOSTCXX="$(which clang++)"
HOSTAR="$(which llvm-ar)"
HOSTRANLIB="$(which llvm-ranlib)"
HOSTAS="$(which llvm-as)"
HOSTSTRIP="$(which llvm-strip)"
HOSTOBJCOPY="$(which llvm-objcopy)"
HOSTOBJDUMP="$(which llvm-objdump)"
HOSTNM="$(which llvm-nm)"

TARGET="x86_64-linux-musl"

function pack() {
    tar -c -C "${SYSROOT}" . | xz -9 -T0 -vvv >"${PROJ_DIR}/sysroot-${1}.tar.xz"
}

function restore() {
    rm -fr "${SYSROOT}"
    mkdir -p "${SYSROOT}"
    env -C "${SYSROOT}" tar -xJf "${PROJ_DIR}/sysroot-${1}.tar.xz"
}
# Idea borrowed from <https://github.com/cbpudding/maplelinux-bootstrap/blob/5484dbc/build-bootstrap.sh#L167>
# Under ISC License
# Copyright (c) 2024-2025 Alexander Hill, Nicholas McDaniel, and Maple Linux Contributors
COMMON_LLVM_CMAKE=(
    "-DCMAKE_ASM_COMPILER_TARGET=${TARGET}"
    "-DCMAKE_BUILD_TYPE=Release"
    "-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON"
    "-DCMAKE_C_COMPILER=${HOSTCC}"
    "-DCMAKE_C_COMPILER_TARGET=${TARGET}"
    "-DCMAKE_CXX_COMPILER=${HOSTCXX}"
    "-DCMAKE_CXX_COMPILER_TARGET=${TARGET}"
    "-DCMAKE_FIND_ROOT_PATH=${SYSROOT}"
    "-DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY"
    "-DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY"
    "-DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=ONLY"
    "-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER"
    "-DCMAKE_INSTALL_PREFIX=${SYSROOT}/usr"
    "-DCMAKE_SYSROOT=${SYSROOT}"
    "-DCMAKE_SYSTEM_NAME=Linux"
    "-DLLVM_HOST_TRIPLE=${TARGET}"
    "-DLLVM_USE_LINKER=lld"
    "-DLLVM_TARGETS_TO_BUILD=X86"
    "-DLLVM_ENABLE_ZSTD=OFF"
    "-DLLVM_ENABLE_ZLIB=OFF"
    "-DLLVM_ENABLE_LIBXML2=OFF"
    "-DCMAKE_CXX_STANDARD=17"
    "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
)