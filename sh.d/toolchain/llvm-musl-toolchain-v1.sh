#!/usr/bin/env bash
set -ue

PROJ_DIR="$(pwd)"
SRC_DIR="${PROJ_DIR}/src"
MUSL_SRC_DIR="${SRC_DIR}/musl-1.2.5"
LLVM_SRC_DIR="${SRC_DIR}/llvm-project-18.1.8.src"
LLVM_CMAKE_SRC_DIR="${LLVM_SRC_DIR}/llvm"
LLVM_COMPILER_RT_SRC_DIR="${LLVM_SRC_DIR}/compiler-rt"

HOSTCC="$(which clang)"
HOSTCXX="$(which clang++)"
HOSTAR="$(which llvm-ar)"
HOSTRANLIB="$(which llvm-ranlib)"
TARGET="x86_64-linux-musl"

rm -fr "${SRC_DIR}"
mkdir -p "${SRC_DIR}"
wget -4 https://musl.libc.org/releases/musl-1.2.5.tar.gz -O "${SRC_DIR}/musl-1.2.5.tar.gz"
env -C "${SRC_DIR}" tar -xf musl-1.2.5.tar.gz


wget -4 https://github.com/llvm/llvm-project/releases/download/llvmorg-18.1.8/llvm-project-18.1.8.src.tar.xz -O "${SRC_DIR}/llvm-project-18.1.8.src.tar.xz"
env -C "${SRC_DIR}" tar -xf llvm-project-18.1.8.src.tar.xz


wget -4 "https://www.kernel.org/pub/linux/kernel/v6.x/linux-6.14.9.tar.xz" -O "${SRC_DIR}/linux-6.14.9.tar.xz"
env -C "${SRC_DIR}" tar -xf linux-6.14.9.tar.xz

rm -fr "${SYSROOT}"
mkdir -p "${SYSROOT}"
env -C "${SYSROOT}" mkdir -pv {dev,proc,run,sys,tmp}
env -C "${SYSROOT}" mkdir -pv usr/{include,lib,libexec,bin}
env -C "${SYSROOT}" ln -sfv usr/bin bin
env -C "${SYSROOT}" ln -sfv usr/bin sbin
env -C "${SYSROOT}" ln -sfv usr/include include
env -C "${SYSROOT}" ln -sfv usr/lib lib
env -C "${SYSROOT}" ln -sfv usr/lib libexec
env -C "${SYSROOT}" ln -sfv usr/lib lib64
env -C "${SYSROOT}" ln -sfv bin usr/sbin
env -C "${SYSROOT}" ln -sfv lib usr/lib64

# Step 1: Build musl headers and Linux kernel headers
env -C \
    "${MUSL_SRC_DIR}" \
    CC="${HOSTCC}" \
    AR="${HOSTAR}" \
    RANLIB="${HOSTRANLIB}" \
    ./configure --prefix="${SYSROOT}"
env -C "${MUSL_SRC_DIR}" make -j"$(nproc)" install-headers

LINUX_SRC_DIR="${SRC_DIR}/linux-6.14.9"
env -C "${LINUX_SRC_DIR}" PATH="${PATH}" make ARCH=x86_64 INSTALL_HDR_PATH="${SYSROOT}/usr" -j"$(nproc)" headers_install

# Back up the current sysroot
tar -cJf "${PROJ_DIR}/sysroot-s1.tar.xz" -C "${SYSROOT}" .

# Step 2: Build clang compiler-rt only, to build musl
LLVM_PHASE1_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-1-clang"
rm -fr "${LLVM_PHASE1_BUILD_DIR}"
mkdir -p "${LLVM_PHASE1_BUILD_DIR}"

cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${SYSROOT}" \
    -DCOMPILER_RT_BAREMETAL_BUILD=ON \
    -DCOMPILER_RT_BUILD_BUILTINS=ON \
    -DCOMPILER_RT_BUILD_LIBFUZZER=OFF \
    -DCOMPILER_RT_BUILD_MEMPROF=OFF \
    -DCOMPILER_RT_BUILD_PROFILE=OFF \
    -DCOMPILER_RT_BUILD_SANITIZERS=OFF \
    -DCOMPILER_RT_BUILD_XRAY=OFF \
    -DCOMPILER_RT_DEFAULT_TARGET_ONLY=ON \
    -DCMAKE_ASM_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_C_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_C_COMPILER="${HOSTCC}" \
    -DCMAKE_CXX_COMPILER="${HOSTCXX}" \
    -DCMAKE_SYSTEM_NAME=Linux \
    -DCMAKE_SYSROOT="${SYSROOT}" \
    -DCMAKE_CXX_COMPILER_WORKS=ON \
    -DCMAKE_C_COMPILER_WORKS=ON \
    -S "${LLVM_COMPILER_RT_SRC_DIR}" \
    -B "${LLVM_PHASE1_BUILD_DIR}"
env -C "${LLVM_PHASE1_BUILD_DIR}" ninja builtins -j"$(nproc)"
env -C "${LLVM_PHASE1_BUILD_DIR}" ninja install-builtins
# Back up the current sysroot
tar -cJf "${PROJ_DIR}/sysroot-s2.tar.xz" -C "${SYSROOT}" .

# Step 3: Build musl with newly built compiler-rt
env -C \
    "${MUSL_SRC_DIR}" \
    CC="${HOSTCC}" \
    AR="${HOSTAR}" \
    RANLIB="${HOSTRANLIB}" \
    LIBCC="${SYSROOT}/lib/linux/libclang_rt.builtins-x86_64.a" \
    ./configure --prefix="${SYSROOT}" \
    CFLAGS="--sysroot=${SYSROOT}"
env -C "${MUSL_SRC_DIR}" make -j"$(nproc)" install
# TODO: Should use links instead of copying
cp "${SYSROOT}/lib/libc.so" "${SYSROOT}/lib/ld-musl-x86_64.so.1"
cp "${SYSROOT}/lib/libc.so" "${SYSROOT}/bin/ldd"
# Back up the musl-installed sysroot
env -C "${SYSROOT}" tar -cJf "${PROJ_DIR}/sysroot-s3.tar.xz" .


# Step 4: Build LLVM RTs: libunwind, libcxxabi, libcxx
LLVM_PHASE4_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-4"
LLVM_RUNTIME_SRC_DIR="${LLVM_SRC_DIR}/runtimes"
rm -fr "${LLVM_PHASE4_BUILD_DIR}"
mkdir -p "${LLVM_PHASE4_BUILD_DIR}"

CFLAGS="-nostdinc++ -rtlib=compiler-rt -Qunused-arguments -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld"
cmake -G "Ninja" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${SYSROOT}/usr" \
    -DCMAKE_C_COMPILER="${HOSTCC}" \
    -DCMAKE_CXX_COMPILER="${HOSTCXX}" \
    -DCMAKE_SYSTEM_NAME=Linux \
    -DLLVM_ENABLE_RUNTIMES="libunwind;libcxxabi;libcxx" \
    -DLIBUNWIND_USE_COMPILER_RT=ON \
    -DLIBUNWIND_ENABLE_SHARED=OFF \
    -DLIBCXX_USE_COMPILER_RT=ON \
    -DLIBCXX_HAS_MUSL_LIBC=ON \
    -DLIBCXX_HAS_ATOMIC_LIB=NO \
    -DLIBCXXABI_USE_LLVM_UNWINDER=ON \
    -DLIBCXXABI_USE_COMPILER_RT=ON \
    -DLIBCXXABI_ENABLE_STATIC_UNWINDER=ON \
    -DCMAKE_SYSROOT="${SYSROOT}" \
    -DLLVM_TARGETS_TO_BUILD=X86 \
    -DCMAKE_C_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_C_COMPILER_WORKS=ON \
    -DCMAKE_CXX_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_CXX_COMPILER_WORKS=ON \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_CXX_STANDARD_REQUIRED=ON \
    -S "${LLVM_RUNTIME_SRC_DIR}" \
    -B "${LLVM_PHASE4_BUILD_DIR}"
cmake --build "${LLVM_PHASE4_BUILD_DIR}" -j"$(nproc)"
cmake --install "${LLVM_PHASE4_BUILD_DIR}"
# Test whether the installed libc++ works
sudo chroot sysroot /bin/ldd /usr/lib/libc++.so.1.0
# Backup the stage 4 sysroot
env -C "${SYSROOT}" tar -cJf "${PROJ_DIR}/sysroot-s4.tar.xz" .

# Step 5: Build compiler-rt with LLVM runtimes linked in
CFLAGS="-nodefaultlibs -nostdinc++ -Qunused-arguments -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld -isystem ${SYSROOT}/usr/include/c++/v1 -stdlib=libc++"
LLVM_PHASE5_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-5"
rm -fr "${LLVM_PHASE5_BUILD_DIR}"
mkdir -p "${LLVM_PHASE5_BUILD_DIR}"

cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${SYSROOT}" \
    -DCOMPILER_RT_BAREMETAL_BUILD=ON \
    -DCOMPILER_RT_BUILD_BUILTINS=ON \
    -DCOMPILER_RT_BUILD_LIBFUZZER=OFF \
    -DCOMPILER_RT_BUILD_MEMPROF=OFF \
    -DCOMPILER_RT_BUILD_PROFILE=OFF \
    -DCOMPILER_RT_BUILD_SANITIZERS=OFF \
    -DCOMPILER_RT_BUILD_XRAY=OFF \
    -DCOMPILER_RT_USE_BUILTINS_LIBRARY=ON \
    -DCOMPILER_RT_DEFAULT_TARGET_ONLY=ON \
    -DCMAKE_ASM_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_C_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_C_COMPILER="${HOSTCC}" \
    -DCMAKE_CXX_COMPILER="${HOSTCXX}" \
    -DCMAKE_C_COMPILER_WORKS=ON \
    -DCMAKE_CXX_FLAGS="${CFLAGS}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_COMPILER_WORKS=ON \
    -DCMAKE_SYSTEM_NAME=Linux \
    -DCMAKE_SYSROOT="${SYSROOT}" \
    -DCMAKE_LIBRARY_PATH="${SYSROOT}/usr/lib" \
    -S "${LLVM_COMPILER_RT_SRC_DIR}" \
    -B "${LLVM_PHASE5_BUILD_DIR}"
env -C "${LLVM_PHASE5_BUILD_DIR}" ninja -j"$(nproc)"
env -C "${LLVM_PHASE5_BUILD_DIR}" ninja install
# Backup the stage 5 sysroot
env -C "${SYSROOT}" tar -cJf "${PROJ_DIR}/sysroot-s5.tar.xz" .

# Step 6: Build clang with LLVM runtimes linked in
CFLAGS="-nodefaultlibs -nostdinc++ -Qunused-arguments -Wl,-lc++ -Wl,-lc++abi -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld -isystem ${SYSROOT}/usr/include/c++/v1 -rtlib=compiler-rt"
LLVM_PHASE6_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-6"
rm -fr "${LLVM_PHASE6_BUILD_DIR}"
mkdir -p "${LLVM_PHASE6_BUILD_DIR}"

cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${SYSROOT}/usr" \
    -DLLVM_ENABLE_LIBCXX=ON \
    -DLLVM_ENABLE_LLD=ON \
    -DCMAKE_SYSTEM_NAME=Linux \
    -DLLVM_ENABLE_ZSTD=OFF \
    -DLLVM_ENABLE_ZLIB=OFF \
    -DCXX_SUPPORTS_CUSTOM_LINKER=true \
    -DLLVM_ENABLE_LIBXML2=OFF \
    -DCMAKE_SYSROOT="${SYSROOT}" \
    -DLLVM_ENABLE_PROJECTS="clang;lld" \
    -DLLVM_TARGETS_TO_BUILD="X86" \
    -DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER \
    -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY \
    -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY \
    -DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=ONLY \
    -DCMAKE_C_COMPILER="${HOSTCC}" \
    -DCMAKE_CXX_COMPILER="${HOSTCXX}" \
    -DCMAKE_C_COMPILER_TARGET="${TARGET}" \
    -DCMAKE_CXX_COMPILER_TARGET="${TARGET}" \
    -DLLVM_HOST_TRIPLE="${TARGET}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_FLAGS="${CFLAGS}" \
    -DCMAKE_SYSROOT="${SYSROOT}" \
    -S "${LLVM_CMAKE_SRC_DIR}" \
    -B "${LLVM_PHASE6_BUILD_DIR}"
env -C "${LLVM_PHASE6_BUILD_DIR}" ninja -j"$(nproc)"
env -C "${LLVM_PHASE6_BUILD_DIR}" ninja install

mkdir -pv "${SYSROOT}/lib/clang/18/lib/linux"
mv "${SYSROOT}/lib/linux" "${SYSROOT}/lib/clang/18/lib/linux"
mv "${SYSROOT}/lib/clang/18/lib/linux/clang_rt.crtbegin-x86_64.o" "${SYSROOT}/lib/clang/18/lib/linux/crtbeginS.o"
mv "${SYSROOT}/lib/clang/18/lib/linux/clang_rt.crtend-x86_64.o" "${SYSROOT}/lib/clang/18/lib/linux/crtendS.o"
mv "${SYSROOT}/lib/clang/18/lib/linux/libclang_rt.builtins.a" "${SYSROOT}/lib/clang/18/lib/linux/libclang_rt.builtins-x86_64.a"

echo "int main() { return 0; }" | sudo chroot "${SYSROOT}" /bin/clang -x c -o a.out - -fuse-ld=lld -rtlib=compiler-rt
sudo chroot "${SYSROOT}" /bin/ldd a.out
printf '#include <stdio.h>\nint main() { printf("Hello, World!\\n"); return 0; }\n' | \
    sudo chroot "${SYSROOT}" /bin/clang -x c -o hello - -fuse-ld=lld -rtlib=compiler-rt -stdlib=libc++
sudo chroot "${SYSROOT}" ./hello

rm -fr "${SYSROOT}"/hello "${SYSROOT}"/a.out
# Backup the step 6 sysroot
env -C "${SYSROOT}" tar -c . | xz -9 -T0 -vvv > "${PROJ_DIR}/sysroot-s6.tar.xz" 

# TODO: Add coreutils, grep, sed, bash, zlib, make, python, cmake, gdb, etc. into the sysroot
# Bash: Use static bash from <https://github.com/robxu9/bash-static/releases/download/5.2.015-1.2.3-2/bash-linux-x86_64>.
