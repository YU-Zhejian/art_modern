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
