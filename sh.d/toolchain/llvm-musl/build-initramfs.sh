#!/usr/bin/env bash
# shellcheck disable=SC1091
set -ue
. vars.sh

make -j"$(nproc)" -C "${LINUX_SRC_DIR}" x86_64_defconfig
make -j"$(nproc)" -C "${LINUX_SRC_DIR}" menuconfig
cp "${LINUX_SRC_DIR}"/.config linux.ini
make -j"$(nproc)" -C "${LINUX_SRC_DIR}" bzImage

BUSYBOX_INSTALL_DIR="${PROJ_DIR}/initramfs"
rm -fr "${BUSYBOX_INSTALL_DIR}"
mkdir -p "${BUSYBOX_INSTALL_DIR}"
env -C "${BUSYBOX_SRC_DIR}" make defconfig
env -C "${BUSYBOX_SRC_DIR}" make menuconfig
cp "${BUSYBOX_SRC_DIR}"/.config busybox.ini

make -j"$(nproc)" -C "${BUSYBOX_SRC_DIR}" \
    CONFIG_PREFIX="${BUSYBOX_INSTALL_DIR}" \
    CC="${HOSTCC}" \
    CFLAGS="--sysroot=${SYSROOT} -isystem ${SYSROOT}/usr/include -rtlib=compiler-rt -Qunused-arguments -fuse-ld=lld -nodefaultlibs -Wl,-lc -static" \
    LD="${HOSTCC}" \
    AS="${HOSTAS}" \
    AR="${HOSTAR}" \
    RANLIB="${HOSTRANLIB}" \
    STRIP="${HOSTSTRIP}" \
    OBJCOPY="${HOSTOBJCOPY}" \
    OBJDUMP="${HOSTOBJDUMP}" \
    NM="${HOSTNM}"
llvm-strip "${BUSYBOX_SRC_DIR}"/busybox_unstripped -o "${BUSYBOX_INSTALL_DIR}"/busybox
mkdir -p "${BUSYBOX_INSTALL_DIR}"/bin
for applet in $(env -C "${BUSYBOX_INSTALL_DIR}" ./busybox --list); do
    env -C "${BUSYBOX_INSTALL_DIR}" ln -sf ../busybox "bin/${applet}"
done

# Pack initramfs
printf '' | cpio -ov -H newc >busybox_initramfs.cpio
for dir in initramfs busybox_initramfs; do
    find "${dir}" |
        sed 's;'"${dir}"';.;' |
        cpio -A -ov -D "$(pwd)/${dir}" -H newc \
            -O busybox_initramfs.cpio
done
gzip -9f busybox_initramfs.cpio

qemu-system-x86_64 \
    -m 2048 \
    -smp 1 \
    -machine pc \
    -kernel "${LINUX_SRC_DIR}"/arch/x86_64/boot/bzImage \
    -initrd busybox_initramfs.cpio.gz \
    -append "${KERNEL_PARAMS}" \
    -nographic \
    -no-reboot
