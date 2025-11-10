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

# MPI run command. Used in MPI-related tests only
# NOTE: Setting this target does NOT enable MPI build!
# Use release-mpi or debug-mpi targets to build with MPI support.
MPIEXEC ?= mpiexec

MPIEXEC_NJOBS ?= 4

# Output directory for build artifacts
OPT_DIR ?= $(CURDIR)/opt

# System architecture
ARCH ?= $(shell uname -m)

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
	mkdir -p $(OPT_DIR)/build_debug
	env -C $(OPT_DIR)/build_debug cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_debug_install/ \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_debug -j$(JOBS)
	cmake --install $(OPT_DIR)/build_debug
	env -C $(OPT_DIR)/build_debug ctest --output-on-failure
	$(OPT_DIR)/build_debug_install/bin/art_modern --help
	$(OPT_DIR)/build_debug_install/bin/art_modern --version

.PHONY: debug-mpi
# debug with MPI
debug-mpi:
	mkdir -p $(OPT_DIR)/build_debug-mpi
	env -C $(OPT_DIR)/build_debug-mpi cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_debug_install-mpi/ \
		-DWITH_MPI=ON \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_debug-mpi -j$(JOBS)
	cmake --install $(OPT_DIR)/build_debug-mpi
	env -C $(OPT_DIR)/build_debug-mpi ctest --output-on-failure
	$(MPIEXEC) -n $(MPIEXEC_NJOBS) $(OPT_DIR)/build_debug_install-mpi/bin/art_modern-mpi --help
	$(MPIEXEC) -n $(MPIEXEC_NJOBS) $(OPT_DIR)/build_debug_install-mpi/bin/art_modern-mpi --version

.PHONY: release
# Generates release build with native optimizations
release:
	mkdir -p $(OPT_DIR)/build_release
	env -C $(OPT_DIR)/build_release cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DCEU_CM_SHOULD_USE_NATIVE=ON \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_release_install/ \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_release -j$(JOBS)
	cmake --install $(OPT_DIR)/build_release
	$(OPT_DIR)/build_release_install/bin/art_modern --help
	$(OPT_DIR)/build_release_install/bin/art_modern --version

.PHONY: release-mpi
# release with MPI
release-mpi:
	mkdir -p $(OPT_DIR)/build_release-mpi
	env -C $(OPT_DIR)/build_release-mpi cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DCEU_CM_SHOULD_USE_NATIVE=ON \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_release_install-mpi/ \
		-DWITH_MPI=ON \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_release-mpi -j$(JOBS)
	cmake --install $(OPT_DIR)/build_release-mpi
	$(MPIEXEC) -n $(MPIEXEC_NJOBS) $(OPT_DIR)/build_release_install-mpi/bin/art_modern-mpi --help
	$(MPIEXEC) -n $(MPIEXEC_NJOBS) $(OPT_DIR)/build_release_install-mpi/bin/art_modern-mpi --version

.PHONY: rel_with_dbg_alpine
# Generates RelWithDebInfo build without native optimizations
# for fully-static-linked build on Alpine Linux
rel_with_dbg_alpine:
	mkdir -p $(OPT_DIR)/build_rel_with_dbg_alpine
	env -C $(OPT_DIR)/build_rel_with_dbg_alpine cmake \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCEU_CM_SHOULD_ENABLE_TEST=OFF \
		-DCEU_CM_SHOULD_USE_NATIVE=OFF \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DBUILD_SHARED_LIBS=OFF \
		-DUSE_MALLOC=NOP \
		-DCMAKE_INSTALL_LIBDIR=bin \
		-DCMAKE_INSTALL_INCLUDEDIR=bin \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_rel_with_dbg_alpine_install/ \
        $(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_rel_with_dbg_alpine -j$(JOBS)
	cmake --install $(OPT_DIR)/build_rel_with_dbg_alpine
	env -C $(OPT_DIR)/build_rel_with_dbg_alpine_install/bin \
		tar cvzf $(OPT_DIR)/build_rel_with_dbg_alpine-$(ARCH).tar.gz \
		art_modern \
		art_profile_builder

.PHONY: fmt
# Run code formatting checks and auto-formatting
# Format the code using [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html), [`sh`](https://github.com/mvdan/sh), [Black](https://black.readthedocs.io/en/stable/), [`cmake-format`](https://cmake-format.readthedocs.io/), and [`dos2unix`](https://www.freebsd.org/cgi/man.cgi?query=dos2unix&sektion=1).
fmt:
	$(BASH) sh.d/fmt.sh

.PHONY: scc
# Run source code counting
# Note that this excludes third-party codes so should be preferred over pure `scc` in project root.
scc:
	$(BASH) sh.d/scc.sh

.PHONY: touch
# Touch all source files to current timestamp
# This **MAY** work when CMake does strange things like compiling the source files again and again.
touch:
	$(BASH) sh.d/touch-all.sh

.PHONY: testsmall
# Run (not necessiarily) small integration tests with debug build. Also tests `art_profile_builder`.
# Env. Flags:
# 
#     - FORMAT_ONLY=1: Stop after testing all output formats is working
#     - NO_FASTQC=1: Do not run FASTQC
testsmall: debug raw_data
	env \
		ART_MODERN_PATH=$(OPT_DIR)/build_debug_install/bin/art_modern \
		APB_PATH=$(OPT_DIR)/build_debug_install/bin/art_profile_builder \
		MPIEXEC="" \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-mpi
# testsmall with MPI
testsmall-mpi: debug-mpi raw_data
	env \
		ART_MODERN_PATH=$(OPT_DIR)/build_debug_install-mpi/bin/art_modern-mpi \
		APB_PATH=$(OPT_DIR)/build_debug_install-mpi/bin/art_profile_builder-mpi \
		MPIEXEC=$(MPIEXEC) \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-release
# Run small tests with release build
testsmall-release: release raw_data
	env \
		ART_MODERN_PATH=$(OPT_DIR)/build_release_install/bin/art_modern \
		APB_PATH=$(OPT_DIR)/build_release_install/bin/art_profile_builder \
		MPIEXEC="" \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-release-mpi
# testsmall-release with MPI
testsmall-release-mpi: release-mpi raw_data
	env \
		ART_MODERN_PATH=$(OPT_DIR)/build_release_install-mpi/bin/art_modern-mpi \
		APB_PATH=$(OPT_DIR)/build_release_install-mpi/bin/art_profile_builder-mpi \
		MPIEXEC=$(MPIEXEC) \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-conda
# Run small tests with conda-installed art_modern
# TODO: Add MPI version when bioconda supports it
testsmall-conda: raw_data
	conda env remove -n _art_modern_bioconda -y || true
	conda create -y -n _art_modern_bioconda -c bioconda -c conda-forge art_modern
	env \
		ART_MODERN_PATH="$(shell conda run -n _art_modern_bioconda type -p art_modern)" \
		APB_PATH="$(shell conda run -n _art_modern_bioconda type -p art_profile_builder)" \
		MPIEXEC="" \
		$(BASH) sh.d/test-small.sh

.PHONY: raw_data
# Download raw data required for tests
raw_data:
	$(MAKE) -C data/raw_data

.PHONY: clean
# Clean all build artifacts
clean:
	rm -fr $(OPT_DIR) tmp build

.PHONY: testbuild
# Test building using diverse conditions
testbuild:
	env MPIEXEC="$(MPIEXEC)" $(PYTHON) sh.d/test-build.py $(CMAKE_FLAGS)

.PHONY: testbuild-mpi
# testbuild with MPI
testbuild-mpi:
	env MPIEXEC="$(MPIEXEC)" $(PYTHON) sh.d/test-build.py --mpi $(CMAKE_FLAGS)

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

.PHONY: serve-doc
# Serve built documentation at <http://localhost:8000>
serve-doc:
	$(PYTHON) -m http.server -d docs/sphinx.d/_build/html

.PHONY: packing
# Create binary packages
packing:
	$(MAKE) -C packing

.PHONY: packing-update-containers
# Update Docker and Singularity containers used for packing
packing-update-containers:
	$(MAKE) -C packing update-containers
