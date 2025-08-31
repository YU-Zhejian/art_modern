CMAKE_FLAGS ?= 
JOBS ?= 40

# build as an alias to debug
.PHONY: build
build: debug

.PHONY: debug
debug:
	mkdir -p opt/build_debug
	env -C opt/build_debug cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		$(CMAKE_FLAGS) \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/build_debug_install/ \
		$(CURDIR)
	cmake --build opt/build_debug -j$(JOBS)
	cmake --install opt/build_debug
	env -C opt/build_debug ctest --output-on-failure
	opt/build_debug_install/bin/art_modern --help
	opt/build_debug_install/bin/art_modern --version
	# cpack --config opt/build_debug/CPackSourceConfig.cmake

.PHONY: release
release:
	mkdir -p opt/build_release
	env -C opt/build_release cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=Release \
		-DCEU_CM_SHOULD_USE_NATIVE=ON \
		$(CMAKE_FLAGS) \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/build_release_install/ \
		$(CURDIR)
	cmake --build opt/build_release -j$(JOBS)
	cmake --install opt/build_release
	opt/build_release_install/bin/art_modern --help
	opt/build_release_install/bin/art_modern --version
	# cpack --config opt/build_release/CPackSourceConfig.cmake

.PHONY: rel_with_dbg_alpine
rel_with_dbg_alpine:
	mkdir -p opt/build_rel_with_dbg_alpine
	env -C opt/build_rel_with_dbg_alpine cmake \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCEU_CM_SHOULD_ENABLE_TEST=OFF \
		-DCEU_CM_SHOULD_USE_NATIVE=OFF \
		-DBUILD_SHARED_LIBS=OFF \
		-DUSE_MALLOC=NOP \
        $(CMAKE_FLAGS) \
		-DCMAKE_INSTALL_LIBDIR=bin \
		-DCMAKE_INSTALL_INCLUDEDIR=bin \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/build_rel_with_dbg_alpine_install/ \
		$(CURDIR)
	cmake --build opt/build_rel_with_dbg_alpine -j$(JOBS)
	cmake --install opt/build_rel_with_dbg_alpine
	env -C opt/build_rel_with_dbg_alpine_install/bin \
		zip -9 -r $(CURDIR)/opt/build_rel_with_dbg_alpine-x86_64.zip \
		art_modern \
		libam_support_lib.a \
		libart_modern_lib.a \
		liblabw_slim_htslib.a \
		libslim_libceu.a \
		libslim_libfmt.a

.PHONY: fmt
fmt:
	bash sh.d/fmt.sh

.PHONY: scc
scc:
	bash sh.d/scc.sh

.PHONY: touch
touch:
	bash sh.d/touch-all.sh

.PHONY: testsmall
testsmall: debug raw_data
	env ART=opt/build_debug_install/bin/art_modern bash sh.d/test_small.sh

.PHONY: testsmall-release
testsmall-release: release raw_data
	env ART=opt/build_release_install/bin/art_modern bash sh.d/test_small.sh

.PHONY: raw_data
raw_data:
	$(MAKE) -C data/raw_data

.PHONY: clean
clean:
	rm -fr opt tmp build

.PHONY: testbuild-child
testbuild-child:
	rm -fr opt/testbuild
	mkdir -p opt/testbuild
	env -C opt/testbuild cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		$(CMAKE_FLAGS) \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/testbuild_install/ \
		$(CURDIR)
	cmake --build opt/testbuild -j$(JOBS)
	cmake --install opt/testbuild
	env -C opt/testbuild ctest --output-on-failure
	opt/testbuild_install/bin/art_modern --help
	opt/testbuild_install/bin/art_modern --version


.PHONY: testbuild
testbuild:
	mkdir -p opt/testbuild
	bash sh.d/test-build.sh

.PHONY: doc
doc:
	$(MAKE) -C docs/sphinx.d
