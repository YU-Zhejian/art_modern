# Makefile for building the art_modern project
# Suffix rules disabled
.SUFFIXES:

# Additional CMake flags for CMake-related tasks
CMAKE_FLAGS ?=

# Number of parallel jobs for building
JOBS ?= $(shell cat /proc/cpuinfo | grep processor | wc -l)

# Shell to use for shell scripts
BASH ?= bash

# Python interpreter to use for Python scripts
PYTHON ?= python3

# Package version, derived from the latest git tag if not set
export PACKAGE_VERSION ?= $(shell git describe --tags --abbrev=0)

$(info Using $(JOBS) parallel jobs for building)
$(info Building for version $(PACKAGE_VERSION))
$(info Using following additional CMake flags: $(CMAKE_FLAGS))

.PHONY: help
help:
	$(PYTHON) $(CURDIR)/sh.d/make2help.py txt < $(CURDIR)/Makefile

.PHONY: build
# build as an alias to debug
build: debug

.PHONY: debug
# Generates debug build with tests enabled
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
# Generates release build with native optimizations
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
# Generates RelWithDebInfo build without native optimizations
# for fully-static-linked build on Alpine Linux
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
# Run code formatting checks and auto-formatting
fmt:
	$(BASH) sh.d/fmt.sh

.PHONY: scc
# Run source code counting
scc:
	$(BASH) sh.d/scc.sh

.PHONY: touch
# Touch all source files to current timestamp
touch:
	$(BASH) sh.d/touch-all.sh

.PHONY: testsmall
# Run small tests with debug build
testsmall: debug raw_data
	env ART=opt/build_debug_install/bin/art_modern $(BASH) sh.d/test-small.sh

.PHONY: testsmall-conda
# Run small tests with conda-installed art_modern
# TODO: This Makefile block requires extensive revision.
testsmall-conda: raw_data
	conda env remove -n _art_modern_bioconda -y || true
	conda create -y -n _art_modern_bioconda -c bioconda -c conda-forge art_modern
	env ART="$(conda run -n _art_modern_bioconda which art_modern)" $(BASH) sh.d/test-small.sh

.PHONY: testsmall-release
# Run small tests with release build
testsmall-release: release raw_data
	env ART=opt/build_release_install/bin/art_modern $(BASH) sh.d/test-small.sh

.PHONY: raw_data
# Download raw data required for tests
raw_data:
	$(MAKE) -C data/raw_data

.PHONY: clean
# Clean all build artifacts
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
# Test building using diverse conditions
testbuild:
	mkdir -p opt/testbuild
	$(BASH) sh.d/test-build.sh

.PHONY: doc
# Build documentation
doc:
	$(PYTHON) $(CURDIR)/sh.d/make2help.py md < $(CURDIR)/Makefile > docs/MakefileTargets.md
	$(MAKE) -C docs/sphinx.d

.PHONY: cleandoc
# Clean and build documentation
cleandoc:
	$(MAKE) -C docs/sphinx.d clean
	$(MAKE) -C docs/sphinx.d

.PHONY: deb
# Build Debian package
#
# **NOTE** This target will fail under WSL
# since the file debian/doc will be regarded as executable under WSL
# which have a different meaning for `debuild`.
deb:
	$(BASH) sh.d/deb.sh
