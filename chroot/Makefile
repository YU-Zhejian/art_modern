
cmake-3.31.4-linux-x86_64.tar.gz:
	wget https://github.com/Kitware/CMake/releases/download/v3.31.4/cmake-3.31.4-linux-x86_64.tar.gz

ubuntu-1804: cmake-3.31.4-linux-x86_64.tar.gz
	sudo debootstrap \
		--variant=buildd \
		--arch=amd64 \
		bionic $(CURDIR)/ubuntu-1804 https://mirrors.tuna.tsinghua.edu.cn/ubuntu/
	sudo cp ubuntu-1804-etc-apt-source.list $(CURDIR)/ubuntu-1804/etc/apt/sources.list
	sudo cp cmake-3.31.4-linux-x86_64.tar.gz ubuntu-1804/opt
	sudo chroot ubuntu-1804 tar -zxvf /opt/cmake-3.31.4-linux-x86_64.tar.gz -C /opt
	sudo chroot ubuntu-1804 apt update
	sudo chroot ubuntu-1804 apt upgrade -y
	sudo chroot ubuntu-1804 apt install -y \
		build-essential \
		g++ \
		binutils \
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
	sudo chroot ubuntu-1804 env -C /root/src/build \
		/opt/cmake-3.31.4-linux-x86_64/bin/cmake \
		/root/src/ \
		-DUSE_THREAD_PARALLEL=NOP \
		-DBOOST_CONFIG_PROVIDED_BY_BOOST=OFF
	sudo chroot ubuntu-1804 env -C /root/src/build make -j40
	sudo chroot ubuntu-1804 env -C /root/src/build \
		/opt/cmake-3.31.4-linux-x86_64/bin/cmake \
		/root/src/ \
		-DUSE_THREAD_PARALLEL=NOP \
		-DBOOST_CONFIG_PROVIDED_BY_BOOST=OFF \
		-DCMAKE_C_COMPILER=clang-5.0 \
		-DCMAKE_CXX_COMPILER=clang++-5.0
	sudo chroot ubuntu-1804 env -C /root/src/build make -j40
