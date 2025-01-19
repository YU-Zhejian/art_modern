

ubuntu-1804:
	sudo debootstrap \
		--variant=buildd \
		--arch=amd64 \
		bionic $(CURDIR)/ubuntu-1804 https://mirrors.tuna.tsinghua.edu.cn/ubuntu/
	sudo cp ubuntu-1804-etc-apt-source.list $(CURDIR)/ubuntu-1804/etc/apt/sources.list
	sudo chroot ubuntu-1804 apt update
	sudo chroot ubuntu-1804 apt upgrade -y
	sudo chroot ubuntu-1804 apt install -y \
		build-essential \
		g++ \
		binutils \
		cmake \
		libboost-all-dev \
		libz-dev \
		make \
		apt-file
	sudo chroot ubuntu-1804 apt-file update
	sudo chroot ubuntu-1804 apt clean

.PHONY: test-build-ubuntu-1804
test-build-ubuntu-1804: ubuntu-1804
	sudo bash update-chroot.sh $(CURDIR)/ubuntu-1804/root/src/
	sudo chroot ubuntu-1804 mkdir -p /root/src/build
	sudo chroot ubuntu-1804 env -C /root/src/build cmake /root/src/ \
		-DUSE_THREAD_PARALLEL=NOP
	sudo chroot ubuntu-1804 env -C /root/src/build make -j40
