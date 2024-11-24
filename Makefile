.PHONY: build
build:
	mkdir -p build
	env -C build cmake -DCMAKE_BUILD_TYPE=Debug -G Ninja ..
	env -C build ninja -j40
	env -C build ctest --output-on-failure
	build/art_modern --help
	mpiexec --verbose -n 5 build/art_modern --version

.PHONY: build_external_htslib
build_external_htslib:
	mkdir -p build_external_htslib
	env -C build_external_htslib cmake \
		-DCMAKE_BUILD_TYPE=Release \
		-DUSE_HTSLIB="hts" \
		-G Ninja ..
	env -C build_external_htslib ninja -j40
	env -C build_external_htslib ctest

.PHONY: fmt
fmt:
	bash fmt.sh

.PHONY: scc
scc:
	bash scc.sh

.PHONY: touch
touch:
	bash touch-all.sh

.PHONY: testsmall
testsmall: build raw_data
	bash test_small.sh

.PHONY: raw_data
raw_data:
	$(MAKE) -C raw_data

.PHONY: build-alpine
build-alpine:
	cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF; ninja

.PHONY: sync-ceu-cm
sync-ceu-cm:
	rm -fr deps/cmake_collections
	cp -rv ../libceu/cmake_collections/ deps/

.PHONY: clean
clean:
	rm -fr build build_profile tmp

.PHONY: profile
profile:
	env -i bash --norc --noprofile profile.sh
