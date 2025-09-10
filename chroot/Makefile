.PHONY: all
all: test-build-ubuntu-1804-gcc test-build-ubuntu-1804-clang alpine-latest-build

cmake-3.17.0-Linux-x86_64.tar.gz:
	wget https://github.com/Kitware/CMake/releases/download/v3.17.0/cmake-3.17.0-Linux-x86_64.tar.gz


debian-oldoldstable: cmake-3.17.0-Linux-x86_64.tar.gz
	sudo debootstrap \
		--variant=buildd \
		--arch=amd64 \
		oldoldstable $(CURDIR)/debian-oldoldstable https://mirrors.ustc.edu.cn/debian/
	sudo cp cmake-3.17.0-Linux-x86_64.tar.gz ubuntu-1804/opt
	sudo cp debian-oldoldstable-etc-apt-source.list $(CURDIR)/debian-oldoldstable/etc/apt/sources.list
	sudo chroot $(CURDIR)/debian-oldoldstable tar -zxvf /opt/cmake-3.17.0-Linux-x86_64.tar.gz -C /opt
	sudo chroot $(CURDIR)/debian-oldoldstable apt update
	sudo chroot $(CURDIR)/debian-oldoldstable apt upgrade -y
	sudo chroot $(CURDIR)/debian-oldoldstable apt install -y \
		build-essential \
		g++ \
		binutils \
		libboost-all-dev \
		libz-dev \
		make \
		apt-file \
		python3 \
		clang++
	sudo chroot $(CURDIR)/debian-oldoldstable apt-file update
	sudo chroot $(CURDIR)/debian-oldoldstable apt clean


ubuntu-1804: cmake-3.17.0-Linux-x86_64.tar.gz
	sudo debootstrap \
		--variant=buildd \
		--arch=amd64 \
		bionic $(CURDIR)/ubuntu-1804 https://www.archive.ubuntu.com/ubuntu/
	sudo cp ubuntu-1804-etc-apt-source.list $(CURDIR)/ubuntu-1804/etc/apt/sources.list
	sudo cp cmake-3.17.0-Linux-x86_64.tar.gz ubuntu-1804/opt
	sudo chroot ubuntu-1804 tar -zxvf /opt/cmake-3.17.0-Linux-x86_64.tar.gz -C /opt
	sudo chroot ubuntu-1804 apt update
	sudo chroot ubuntu-1804 apt upgrade -y
	sudo chroot ubuntu-1804 apt install -y \
		build-essential \
		g++ \
		binutils \
		libboost-all-dev \
		libz-dev \
		make \
		apt-file \
		python3 \
		clang++-5
	sudo chroot ubuntu-1804 apt-file update
	sudo chroot ubuntu-1804 apt clean

.PHONY: test-build-ubuntu-1804-gcc
test-build-ubuntu-1804-gcc: ubuntu-1804
	sudo bash update-chroot.sh $(CURDIR)/ubuntu-1804/root/src/
	sudo chroot ubuntu-1804 rm -fr /root/src/build
	sudo chroot ubuntu-1804 mkdir -p /root/src/build
	sudo chroot ubuntu-1804 env -C /root/src/build \
		/opt/cmake-3.17.0-Linux-x86_64/bin/cmake \
		/root/src/ \
		-DUSE_THREAD_PARALLEL=BS \
		-DBOOST_CONFIG_PROVIDED_BY_BOOST=OFF \
        -DCMAKE_C_COMPILER=gcc \
        -DCMAKE_CXX_COMPILER=g++
	sudo chroot ubuntu-1804 env -C /root/src/build make -j40
	sudo chroot ubuntu-1804 env -C /root/src/build \
		/opt/cmake-3.17.0-Linux-x86_64/bin/ctest --output-on-failure
	sudo chroot ubuntu-1804 /root/src/build/art_modern --help
	sudo chroot ubuntu-1804 /root/src/build/art_modern --version

.PHONY: test-build-ubuntu-1804-clang
test-build-ubuntu-1804-clang: ubuntu-1804
	sudo bash update-chroot.sh $(CURDIR)/ubuntu-1804/root/src/
	sudo chroot ubuntu-1804 rm -fr /root/src/build
	sudo chroot ubuntu-1804 mkdir -p /root/src/build
	sudo chroot ubuntu-1804 env -C /root/src/build \
		/opt/cmake-3.17.0-Linux-x86_64/bin/cmake \
		/root/src/ \
		-DUSE_THREAD_PARALLEL=BS \
		-DBOOST_CONFIG_PROVIDED_BY_BOOST=OFF \
		-DCMAKE_C_COMPILER=clang-5.0 \
		-DCMAKE_CXX_COMPILER=clang++-5.0
	sudo chroot ubuntu-1804 env -C /root/src/build make -j40
	sudo chroot ubuntu-1804 env -C /root/src/build \
		/opt/cmake-3.17.0-Linux-x86_64/bin/ctest --output-on-failure
	sudo chroot ubuntu-1804 /root/src/build/art_modern --help
	sudo chroot ubuntu-1804 /root/src/build/art_modern --version

alpine-latest:
	git clone https://gitlab.alpinelinux.org/alpine/apk-tools apk-tools
	env -C apk-tools meson setup -Dc_link_args="-static" -Dprefer_static=true -Ddefault_library=static build
	env -C apk-tools ninja -C build src/apk
	sudo $(CURDIR)/apk-tools/build/src/apk \
		-X https://mirrors.tuna.tsinghua.edu.cn/alpine/latest-stable/main \
		-U --allow-untrusted \
		-p alpine-latest \
		--initdb add alpine-base
	sudo cp /etc/resolv.conf alpine-latest/etc/resolv.conf
	sudo cp alpine-latest-etc-apk-repositories alpine-latest/etc/apk/repositories
	sudo mount -o bind /dev alpine-latest/dev
	sudo mount -t proc none alpine-latest/proc
	sudo mount -o bind /sys alpine-latest/sys
	sudo chroot alpine-latest apk update
	sudo chroot alpine-latest apk add \
        g++ \
        binutils \
        boost-dev \
        zlib-dev \
        make \
        python3 \
        cmake \
        zlib-static \
		icu-static \
		coreutils \
		boost1.84-static
	# FIXME: The boost version may subject to change.

.PHONY: alpine-latest-build
alpine-latest-build: alpine-latest
	sudo bash update-chroot.sh $(CURDIR)/alpine-latest/root/src/
	sudo chroot $(CURDIR)/alpine-latest/ rm -fr /root/src/build
	sudo chroot $(CURDIR)/alpine-latest/ mkdir -p /root/src/build
	sudo chroot $(CURDIR)/alpine-latest/ \
		env -C /root/src/build cmake \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCEU_CM_SHOULD_ENABLE_TEST=OFF \
		-DCEU_CM_SHOULD_USE_NATIVE=OFF \
		-DBUILD_SHARED_LIBS=OFF \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
		/root/src
	sudo chroot $(CURDIR)/alpine-latest/ \
		cmake --build /root/src/build -j40
	sudo cp $(CURDIR)/alpine-latest/root/src/build/art_modern $(CURDIR)/art_modern-alpine
	sudo chown $(shell id -u):$(shell id -g) $(CURDIR)/art_modern-alpine
