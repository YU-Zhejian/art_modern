#!/usr/bin/env bash
# shellcheck disable=SC1091
set -ue

. vars.sh

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

env -C "${LINUX_SRC_DIR}" \
    make \
    -j"$(nproc)" \
    headers_install \
    ARCH=x86_64 \
    INSTALL_HDR_PATH="${SYSROOT}/usr"

# Back up the current sysroot
pack s1

# Step 2: Build clang compiler-rt only, to build musl
LLVM_PHASE1_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-1"
rm -fr "${LLVM_PHASE1_BUILD_DIR}"
mkdir -p "${LLVM_PHASE1_BUILD_DIR}"

cmake -G Ninja "${COMMON_LLVM_CMAKE[@]}" \
    -DCOMPILER_RT_BAREMETAL_BUILD=ON \
    -DCOMPILER_RT_BUILD_BUILTINS=ON \
    -DCOMPILER_RT_BUILD_LIBFUZZER=OFF \
    -DCOMPILER_RT_BUILD_MEMPROF=OFF \
    -DCOMPILER_RT_BUILD_PROFILE=OFF \
    -DCOMPILER_RT_BUILD_SANITIZERS=OFF \
    -DCOMPILER_RT_BUILD_XRAY=OFF \
    -DCOMPILER_RT_DEFAULT_TARGET_ONLY=ON \
    -DCMAKE_CXX_COMPILER_WORKS=ON \
    -DCMAKE_C_COMPILER_WORKS=ON \
    -S "${LLVM_COMPILER_RT_SRC_DIR}" \
    -B "${LLVM_PHASE1_BUILD_DIR}"
env -C "${LLVM_PHASE1_BUILD_DIR}" ninja builtins -j"$(nproc)"
env -C "${LLVM_PHASE1_BUILD_DIR}" ninja install-builtins
# Back up the current sysroot
pack s2

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
env -C "${SYSROOT}/lib/" ln -sfv libc.so ld-musl-x86_64.so.1
env -C "${SYSROOT}/bin/" ln -sfv ../lib/libc.so ldd
# Back up the musl-installed sysroot
pack s3

# Step 4: Build LLVM RTs: libunwind, libcxxabi, libcxx
LLVM_PHASE4_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-4"
LLVM_RUNTIME_SRC_DIR="${LLVM_SRC_DIR}/runtimes"
rm -fr "${LLVM_PHASE4_BUILD_DIR}"
mkdir -p "${LLVM_PHASE4_BUILD_DIR}"

CFLAGS="-nostdinc++ -rtlib=compiler-rt -Qunused-arguments -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld -w"
cmake -G "Ninja" "${COMMON_LLVM_CMAKE[@]}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_C_COMPILER_WORKS=ON \
    -DCMAKE_CXX_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_COMPILER_WORKS=ON \
    -DLLVM_ENABLE_RUNTIMES="libunwind;libcxxabi;libcxx" \
    -DLIBUNWIND_USE_COMPILER_RT=ON \
    -DLIBUNWIND_ENABLE_SHARED=OFF \
    -DLIBCXX_USE_COMPILER_RT=ON \
    -DLIBCXX_HAS_MUSL_LIBC=ON \
    -DLIBCXX_HAS_ATOMIC_LIB=NO \
    -DLIBCXXABI_USE_LLVM_UNWINDER=ON \
    -DLIBCXXABI_USE_COMPILER_RT=ON \
    -S "${LLVM_RUNTIME_SRC_DIR}" \
    -B "${LLVM_PHASE4_BUILD_DIR}"
cmake --build "${LLVM_PHASE4_BUILD_DIR}" -j"$(nproc)"
cmake --install "${LLVM_PHASE4_BUILD_DIR}"
# Test whether the installed libc++ works
sudo chroot --userspec="$(id -u):$(id -g)" sysroot /bin/ldd /usr/lib/libc++.so.1.0
# Backup the stage 4 sysroot
pack s4

# Step 5: Build compiler-rt with LLVM runtimes linked in

# Substep: Build and install libexecinfo
CFLAGS="-nodefaultlibs -Qunused-arguments -isystem ${SYSROOT}/usr/include/ -w --sysroot=${SYSROOT} --target=${TARGET}"
LDFLAGS="-Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld"
env -C "${LIBEXECINFO_SRC_DIR}" \
    make install -j"$(nproc)" CC="${HOSTCC}" \
    CFLAGS="${CFLAGS}" \
    LDFLAGS="${LDFLAGS}" \
    PREFIX="${SYSROOT}/usr"
sudo chroot --userspec="$(id -u):$(id -g)" sysroot /bin/ldd /usr/lib/libexecinfo.so.1

CFLAGS="-nodefaultlibs -nostdinc++ -Qunused-arguments -Wl,-lexecinfo -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld -isystem ${SYSROOT}/usr/include/c++/v1 -isystem ${SYSROOT}/usr/include/ -stdlib=libc++ -w"
LLVM_PHASE5_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-5"
rm -fr "${LLVM_PHASE5_BUILD_DIR}"
mkdir -p "${LLVM_PHASE5_BUILD_DIR}"

cmake -G Ninja "${COMMON_LLVM_CMAKE[@]}" \
    -DCOMPILER_RT_BAREMETAL_BUILD=ON \
    -DCOMPILER_RT_BUILD_BUILTINS=ON \
    -DCOMPILER_RT_BUILD_LIBFUZZER=OFF \
    -DCOMPILER_RT_BUILD_MEMPROF=ON \
    -DCOMPILER_RT_BUILD_PROFILE=ON \
    -DCOMPILER_RT_BUILD_SANITIZERS=ON \
    -DCOMPILER_RT_BUILD_XRAY=ON \
    -DCOMPILER_RT_USE_BUILTINS_LIBRARY=ON \
    -DCOMPILER_RT_DEFAULT_TARGET_ONLY=ON \
    -DCMAKE_C_COMPILER_WORKS=ON \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_COMPILER_WORKS=ON \
    -DCMAKE_CXX_FLAGS="${CFLAGS}" \
    -DCMAKE_LIBRARY_PATH="${SYSROOT}/usr/lib" \
    -S "${LLVM_COMPILER_RT_SRC_DIR}" \
    -B "${LLVM_PHASE5_BUILD_DIR}"
env -C "${LLVM_PHASE5_BUILD_DIR}" ninja -j"$(nproc)"
env -C "${LLVM_PHASE5_BUILD_DIR}" ninja install
# Backup the stage 5 sysroot
pack s5

# Step 6: Build clang with LLVM runtimes linked in
CFLAGS="-nodefaultlibs -nostdinc++ -Qunused-arguments -Wl,-lc++ -Wl,-lc++abi -Wl,-lexecinfo -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld -isystem ${SYSROOT}/usr/include/c++/v1 -isystem ${SYSROOT}/usr/include/ -rtlib=compiler-rt -unwindlib=libunwind"
LLVM_PHASE6_BUILD_DIR="${PROJ_DIR}/build/llvm-phase-6"
rm -fr "${LLVM_PHASE6_BUILD_DIR}"
mkdir -p "${LLVM_PHASE6_BUILD_DIR}"

cmake -G Ninja "${COMMON_LLVM_CMAKE[@]}" \
    -DLLVM_ENABLE_LIBCXX=ON \
    -DCXX_SUPPORTS_CUSTOM_LINKER=true \
    -DLLVM_ENABLE_PROJECTS="clang;lld" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_FLAGS="${CFLAGS}" \
    -DCLANG_DEFAULT_CXX_STDLIB=libc++ \
    -DCLANG_DEFAULT_OBJCOPY=llvm-objcopy \
    -DCLANG_DEFAULT_LINKER=lld \
    -DCLANG_DEFAULT_RTLIB=compiler-rt \
    -DCLANG_DEFAULT_UNWINDLIB=libunwind \
    -S "${LLVM_CMAKE_SRC_DIR}" \
    -B "${LLVM_PHASE6_BUILD_DIR}"
env -C "${LLVM_PHASE6_BUILD_DIR}" ninja -j"$(nproc)"
env -C "${LLVM_PHASE6_BUILD_DIR}" ninja install

mkdir -pv "${SYSROOT}/lib/clang/18/lib"
mv -v "${SYSROOT}/lib/linux" "${SYSROOT}/lib/clang/18/lib/"
mv -v "${SYSROOT}/lib/clang/18/lib/linux/clang_rt.crtbegin-x86_64.o" "${SYSROOT}/lib/clang/18/lib/linux/crtbeginS.o"
mv -v "${SYSROOT}/lib/clang/18/lib/linux/clang_rt.crtend-x86_64.o" "${SYSROOT}/lib/clang/18/lib/linux/crtendS.o"

# Test the final toolchain
echo "int main() { return 0; }" | sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /bin/clang -x c -o a.out -
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /bin/ldd a.out
printf '#include <stdio.h>\nint main() { printf("Hello, World!\\n"); return 0; }\n' |
    sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /bin/clang -x c -o hello -
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" ./hello
rm -fr "${SYSROOT}"/hello "${SYSROOT}"/a.out
# Backup the step 6 sysroot
pack s6

# Step 7: Build GNU Coreutils
env -C "${COREUTILS_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --with-openssl=no \
    --with-linux-crypto \
    --without-selinux \
    --without-libgmp \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld" \
    CPPFLAGS="--sysroot=${SYSROOT} -isystem ${SYSROOT}/usr/include" \
    CC="${HOSTCC}"
env -C "${COREUTILS_SRC_DIR}" make -j"$(nproc)" install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/ldd /usr/bin/ls
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/ls -lFh
pack s7

# Step 8: Build Bash
env -C "${BASH_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --without-bash-malloc \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld" \
    CPPFLAGS="--sysroot=${SYSROOT} -isystem ${SYSROOT}/usr/include" \
    CC="${HOSTCC}"
env -C "${BASH_SRC_DIR}" make -j"$(nproc)" install
cp -v "${SYSROOT}/usr/bin/bash" "${SYSROOT}/usr/bin/sh"
sudo env -i "$(which chroot)" --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/bash -li
pack s8

# Step 9: Build util-linux
env -C "${UTIL_LINUX_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --disable-nls \
    --disable-pylibmount \
    --disable-liblastlog2 \
    --disable-rfkill \
    --disable-bash-completion \
    --disable-makeinstall-chown \
    --disable-makeinstall-setuid \
    --disable-makeinstall-tty-setgid \
    --without-python \
    --without-util \
    --without-tinfo \
    --without-udev \
    --without-ncursesw \
    --without-ncurses \
    --without-btrfs \
    --without-systemd \
    --disable-asciidoc \
    --disable-poman \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld -Wl,-nodefaultlibs" \
    CC="${HOSTCC}" \
    CXX="${HOSTCXX}" \
    CXXFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc++ -Wl,-lc++abi -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -nostdinc++ -isystem ${SYSROOT}/usr/include/c++/v1 -isystem ${SYSROOT}/usr/include/ -stdlib=libc++ -w -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld"
env -C "${UTIL_LINUX_SRC_DIR}" \
    make -j"$(nproc)" AM_DEFAULT_VERBOSITY=1 \
    LIBS="-XCClinker -nodefaultlibs" \
    install
sudo env -i "$(which chroot)" --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/mount --help
pack s9

# Step 10: Build grep
# AWK not build since not required in the minimal toolchain
env -C "${GREP_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --disable-threads \
    --disable-nls \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld" \
    CC="${HOSTCC}"
env -C "${GREP_SRC_DIR}" make -j"$(nproc)"
env -C "${GREP_SRC_DIR}" make install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/grep --version

env -C "${SED_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --disable-threads \
    --disable-nls \
    --without-selinux \
    --disable-i18n \
    --disable-acl \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld" \
    CC="${HOSTCC}"
env -C "${SED_SRC_DIR}" make -j"$(nproc)"
env -C "${SED_SRC_DIR}" make install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/sed --version
pack s10

# Stage 11: Add zlib, xz, make, pkgconf
env -C "${MAKE_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --disable-nls \
    --without-guile \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld  -D_XOPEN_SOURCE=500" \
    CC="${HOSTCC}"
env -C "${MAKE_SRC_DIR}" make -j"$(nproc)" install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/make --version

env -C "${ZLIB_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr"
env -C "${ZLIB_SRC_DIR}" make -j"$(nproc)" \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind" \
    CC="${HOSTCC}" \
    install

env -C "${PKGCONF_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99" \
    CC="${HOSTCC}"
env -C "${PKGCONF_SRC_DIR}" make -j"$(nproc)" \
    AM_DEFAULT_VERBOSITY=1 \
    LIBS='-XCClinker -nodefaultlibs -XCClinker -fuse-ld=lld -XCClinker -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1' \
    install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" \
    /usr/bin/env PKG_CONFIG_LIBDIR=/usr/lib/pkgconfig \
    /usr/bin/pkgconf --list-all

env -C "${XZ_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --disable-nls \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld" \
    CC="${HOSTCC}"
env -C "${XZ_SRC_DIR}" make -j"$(nproc)" \
    LIBS='-XCClinker -nodefaultlibs -XCClinker -fuse-ld=lld -XCClinker -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1' install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/xz --version

env -C "${LIBFFI_SRC_DIR}" ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --with-sysroot="${SYSROOT}" \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -Wc,-nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld" \
    CC="${HOSTCC}" \
    CXX="${HOSTCXX}" \
    CXXFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -Wc,-nodefaultlibs -Wl,-lc++ -Wl,-lc++abi -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -nostdinc++ -isystem ${SYSROOT}/usr/include/c++/v1 -isystem ${SYSROOT}/usr/include/ -stdlib=libc++ -w -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld"

env -C "${LIBFFI_SRC_DIR}" make -j"$(nproc)"
env -C "${LIBFFI_SRC_DIR}" make install

env -C "${BZIP2_SRC_DIR}" make \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind" \
    CC="${HOSTCC}" \
    PREFIX="${SYSROOT}/usr" \
    install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/bzip2 --version
pack s11

# Step 12: Build Python
env -C "${PYTHON_SRC_DIR}" CONFIG_SITE="${PROJ_DIR}"/py.config.site \
    ./configure \
    --prefix="${SYSROOT}/usr" \
    --host="${TARGET}" \
    --build="x86_64-pc-linux-gnu" \
    --disable-ipv6 \
    --disable-static \
    --enable-big-digits \
    --without-pydebug \
    --without-pymalloc \
    --without-doc-strings \
    --without-c-locale-coercion \
    --without-readline \
    --without-ensurepip \
    --with-build-python="${HOST_PYTHON}" \
    CC="${HOSTCC}" \
    CFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -w -std=gnu99 -rtlib=compiler-rt -unwindlib=libunwind -fuse-ld=lld -fPIC" \
    CPPFLAGS="--sysroot=${SYSROOT} -isystem ${SYSROOT}/usr/include" \
    LDFLAGS="--sysroot=${SYSROOT} -Qunused-arguments -nodefaultlibs -Wl,-lc -Wl,--dynamic-linker=/lib/ld-musl-x86_64.so.1 -fuse-ld=lld"
env -C "${PYTHON_SRC_DIR}" make -j"$(nproc)" install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/python3 --version
pack s12

# Step 13: Build CMake
mkdir -p "${SYSROOT}/usr/src"
cp -r "${CMAKE_SRC_DIR}" "${SYSROOT}/usr/src/cmake"
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" \
    env -C /usr/src/cmake PATH="/usr/bin" MAKE="/usr/bin/make" \
    bash ./bootstrap \
    --no-system-libs \
    --prefix="/usr" \
    --no-qt-gui \
    --parallel="$(nproc)" \
    -- -DCMAKE_USE_OPENSSL=OFF
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" \
    env -C /usr/src/cmake PATH="/usr/bin" MAKE="/usr/bin/make" \
    /usr/bin/make -j"$(nproc)"
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" \
    env -C /usr/src/cmake PATH="/usr/bin" MAKE="/usr/bin/make" \
    /usr/bin/make install
sudo chroot --userspec="$(id -u):$(id -g)" "${SYSROOT}" /usr/bin/cmake --version
pack s13

# Step 14: TODO: Build Boost

mksquashfs "${SYSROOT}" rootfs.squashfs -comp zstd -noappend
# Create data qcow2 image
qemu-img create -f qcow2 data.qcow2 40G

qemu-system-x86_64 \
    -m 2048 \
    -smp 1 \
    -machine pc \
    -kernel "${LINUX_SRC_DIR}"/arch/x86_64/boot/bzImage \
    -initrd busybox_initramfs.cpio.gz \
    -drive file=rootfs.squashfs,format=raw,if=ide,index=0 \
    -drive file=data.qcow2,format=qcow2,if=ide,index=1 \
    -append "${KERNEL_PARAMS}" \
    -nographic \
    -no-reboot
