CMAKE_FLAGS ?= 
JOBS ?= $(shell cat /proc/cpuinfo | grep processor | wc -l)
BASH ?= bash

export PACKAGE_VERSION ?= $(shell git describe --tags --abbrev=0)

$(info Using $(JOBS) parallel jobs for building)
$(info Building for version $(PACKAGE_VERSION))
$(info Using following additional CMake flags: $(CMAKE_FLAGS))

.PHONY: help
help:
	@echo "debug release rel_with_dbg_alpine fmt scc touch testsmall testsmall-conda testsmall-release raw_data clean testbuild doc cleandoc deb"

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
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/build_debug_install/ \
		$(CMAKE_FLAGS) \
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
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/build_release_install/ \
		$(CMAKE_FLAGS) \
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
		-DCMAKE_INSTALL_LIBDIR=bin \
		-DCMAKE_INSTALL_INCLUDEDIR=bin \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)/opt/build_rel_with_dbg_alpine_install/ \
        $(CMAKE_FLAGS) \
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
	$(BASH) sh.d/fmt.sh

.PHONY: scc
scc:
	$(BASH) sh.d/scc.sh

.PHONY: touch
touch:
	$(BASH) sh.d/touch-all.sh

.PHONY: testsmall
testsmall: debug raw_data
	env ART=opt/build_debug_install/bin/art_modern $(BASH) sh.d/test_small.sh

.PHONY: testsmall-conda
testsmall-conda: raw_data
	# TODO: This Makefile block requires extensive revision.
	conda env remove -n _art_modern_bioconda -y || true
	conda create -y -n _art_modern_bioconda -c bioconda -c conda-forge art_modern
	env ART="$(conda run -n _art_modern_bioconda which art_modern)" $(BASH) sh.d/test_small.sh

.PHONY: testsmall-release
testsmall-release: release raw_data
	env ART=opt/build_release_install/bin/art_modern $(BASH) sh.d/test_small.sh

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
	$(BASH) sh.d/test-build.sh

.PHONY: doc
doc:
	$(MAKE) -C docs/sphinx.d

.PHONY: cleandoc
cleandoc:
	$(MAKE) -C docs/sphinx.d clean
	$(MAKE) -C docs/sphinx.d

.PHONY: deb
deb:
	# NOTE: This target will fail under WSL
	# since the file debian/doc will be regarded as executable under WSL
	# which have a different meaning for debuild.
	$(BASH) sh.d/deb.sh
