# Makefile for building the art_modern project
# Suffix rules disabled
.SUFFIXES:

# Add -R to make to clear buitin variables
MAKEFLAGS += -R

# Additional CMake flags for CMake-related tasks
CMAKE_FLAGS ?=

# Number of parallel jobs for building.
# Defaults to number of CPU cores detected using `nproc`.
# Guarded ':=' not '?=': '?=' is recursive, re-running `nproc` on every reference.
ifeq ($(origin JOBS),undefined)
JOBS := $(shell nproc)
endif

# Shell to use for shell scripts
BASH ?= bash

# Python interpreter to use for Python scripts.
ifeq ($(origin PIXI_PYTHON),undefined)
PIXI_PYTHON := $(shell pixi run --frozen which python)
endif

# Prefix for recipes needing the interpreter: an empty PIXI_PYTHON would otherwise leave them
# running the script as a bare command.
REQUIRE_PIXI_PYTHON = $(if $(strip $(PIXI_PYTHON)),,$(error PIXI_PYTHON is empty: 'pixi run --frozen which python' produced no output. Is pixi installed, and is the default environment solved?))

# MPI launcher command. Used in MPI-related tests only
# NOTE: Setting this target does NOT enable MPI build!
# Use release-mpi or debug-mpi targets to build with MPI support.
MPIEXEC ?= mpiexec

# Output directory for build artifacts
export OPT_DIR ?= $(CURDIR)/opt

# System architecture
ifeq ($(origin ARCH),undefined)
ARCH := $(shell uname -m)
endif

# Package version.
ifeq ($(origin PACKAGE_VERSION),undefined)
PACKAGE_VERSION := $(shell tr -d '[:space:]' < version.txt)
endif
export PACKAGE_VERSION

$(info Using $(JOBS) parallel jobs for building)
$(info Building for version $(PACKAGE_VERSION))
$(info Using following additional CMake flags: $(CMAKE_FLAGS))

# Abort if we're inside a PIXI environment, as this Makefile is not intended to be used inside PIXI.
ifneq ($(PIXI_ENVIRONMENT_NAME),)
$(error This Makefile is not intended to be used inside PIXI. Please run it outside of PIXI environment.)
endif

.PHONY: help
# Generate help message from Makefile comments
help:
	$(REQUIRE_PIXI_PYTHON)$(PIXI_PYTHON) $(CURDIR)/sh.d/make2help.py txt < $(CURDIR)/Makefile

.PHONY: dumpvar
# List all variables using .VARIABLES
dumpvar:
	@$(foreach v,$(.VARIABLES),$(info $(v)=$($(v))))


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
	env -C $(OPT_DIR)/build_debug ctest --output-on-failure
	cmake --install $(OPT_DIR)/build_debug

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
	env -C $(OPT_DIR)/build_debug-mpi ctest --output-on-failure
	cmake --install $(OPT_DIR)/build_debug-mpi

.PHONY: release
# Generates release build with native optimizations
release:
	mkdir -p $(OPT_DIR)/build_release
	env -C $(OPT_DIR)/build_release cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DCEU_CM_SHOULD_USE_NATIVE=ON \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		-DBUILD_ART_MODERN_BENCHMARKS=ON \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_release_install/ \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_release -j$(JOBS)
	env -C $(OPT_DIR)/build_release ctest --output-on-failure
	cmake --install $(OPT_DIR)/build_release

.PHONY: release-mpi
# release with MPI
release-mpi:
	mkdir -p $(OPT_DIR)/build_release-mpi
	env -C $(OPT_DIR)/build_release-mpi cmake \
		-Wdev -Wdeprecated --warn-uninitialized \
		-DCMAKE_BUILD_TYPE=RelWithDebInfo \
		-DCMAKE_VERBOSE_MAKEFILE=ON \
		-DCEU_CM_SHOULD_USE_NATIVE=ON \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		-DBUILD_ART_MODERN_BENCHMARKS=ON \
		-DCMAKE_INSTALL_LIBDIR=lib/art_modern/lib \
		-DCMAKE_INSTALL_INCLUDEDIR=include/art_modern/include \
		-DCMAKE_INSTALL_PREFIX=$(OPT_DIR)/build_release_install-mpi/ \
		-DWITH_MPI=ON \
		$(CMAKE_FLAGS) \
		$(CURDIR)
	cmake --build $(OPT_DIR)/build_release-mpi -j$(JOBS)
	env -C $(OPT_DIR)/build_release-mpi ctest --output-on-failure
	cmake --install $(OPT_DIR)/build_release-mpi

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
		art_profile_builder \
		am_compress

.PHONY: fmt
# Run code formatting checks and auto-formatting
# Format the code using [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html), [`sh`](https://github.com/mvdan/sh), [Black](https://black.readthedocs.io/en/stable/), [`cmake-format`](https://cmake-format.readthedocs.io/), and [`dos2unix`](https://www.freebsd.org/cgi/man.cgi?query=dos2unix&sektion=1).
fmt:
	pixi run --frozen -e fmt $(BASH) sh.d/fmt.sh

.PHONY: lint
# Run static analysis using [`clang-tidy`](https://clang.llvm.org/extra/clang-tidy/), [Pyright](https://microsoft.github.io/pyright/) and [ShellCheck](https://www.shellcheck.net/).
# Configures `cmake-build-clangd` for the compilation database, which is also what `clangd` reads.
lint:
	cmake -S $(CURDIR) -B $(CURDIR)/cmake-build-clangd \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
		-DCEU_CM_SHOULD_ENABLE_TEST=ON \
		$(CMAKE_FLAGS)
	pixi run --frozen -e lint $(BASH) sh.d/lint.sh

.PHONY: scc
# Run source code counting
# Note that this excludes third-party codes so should be preferred over pure `scc` in project root.
scc:
	$(BASH) sh.d/scc.sh

.PHONY: testsmall
# Alias to testsmall-debug
testsmall: testsmall-debug

.PHONY: testsmall-mpi
# Alias to testsmall-debug-mpi
testsmall-mpi: testsmall-debug-mpi

.PHONY: testsmall-debug
# Run (not necessiarily) small integration tests with debug build. Also tests `art_profile_builder`.
# Env. Flags:
# 
#     - FORMAT_ONLY=1: Stop after testing all output formats is working
#     - NO_FASTQC=1: Do not run FASTQC
testsmall-debug: debug raw_data
	pixi run --frozen env \
		ART_MODERN_PATH=$(OPT_DIR)/build_debug_install/bin/art_modern \
		APB_PATH=$(OPT_DIR)/build_debug_install/bin/art_profile_builder \
		AMC_PATH=$(OPT_DIR)/build_debug_install/bin/am_compress \
		WITH_MPI=0 \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-debug-mpi
# testsmall with MPI
testsmall-debug-mpi: debug-mpi raw_data
	pixi run --frozen env \
		ART_MODERN_PATH=$(OPT_DIR)/build_debug_install-mpi/bin/art_modern-mpi \
		APB_PATH=$(OPT_DIR)/build_debug_install-mpi/bin/art_profile_builder-mpi \
		AMC_PATH=$(OPT_DIR)/build_debug_install-mpi/bin/am_compress-mpi \
		WITH_MPI=1 \
		MPIEXEC=$(MPIEXEC) \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-release
# Run small tests with release build
testsmall-release: release raw_data
	pixi run --frozen env \
		ART_MODERN_PATH=$(OPT_DIR)/build_release_install/bin/art_modern \
		APB_PATH=$(OPT_DIR)/build_release_install/bin/art_profile_builder \
		AMC_PATH=$(OPT_DIR)/build_release_install/bin/am_compress \
		WITH_MPI=0 \
		$(BASH) sh.d/test-small.sh

.PHONY: testsmall-release-mpi
# testsmall-release with MPI
testsmall-release-mpi: release-mpi raw_data
	pixi run --frozen env \
		ART_MODERN_PATH=$(OPT_DIR)/build_release_install-mpi/bin/art_modern-mpi \
		APB_PATH=$(OPT_DIR)/build_release_install-mpi/bin/art_profile_builder-mpi \
		AMC_PATH=$(OPT_DIR)/build_release_install-mpi/bin/am_compress-mpi \
		WITH_MPI=1 \
		MPIEXEC=$(MPIEXEC) \
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
	$(REQUIRE_PIXI_PYTHON)$(PIXI_PYTHON) sh.d/test-build.py $(CMAKE_FLAGS)

.PHONY: testbuild-mpi
# testbuild with MPI
testbuild-mpi:
	$(REQUIRE_PIXI_PYTHON)env MPIEXEC="$(MPIEXEC)" $(PIXI_PYTHON) sh.d/test-build.py --mpi $(CMAKE_FLAGS)

.PHONY: testbuild-small
# Test building using diverse conditions with testsmall configurations
testbuild-small:
	$(REQUIRE_PIXI_PYTHON)$(PIXI_PYTHON) sh.d/test-build.py --small $(CMAKE_FLAGS)

.PHONY: testbuild-small-mpi
# testbuild with MPI with testsmall configurations
testbuild-small-mpi:
	$(REQUIRE_PIXI_PYTHON)env MPIEXEC="$(MPIEXEC)" $(PIXI_PYTHON) sh.d/test-build.py --small --mpi $(CMAKE_FLAGS)

.PHONY: doc
# Build documentation
# Deliberately not $(PIXI_PYTHON): Sphinx below needs the `doc` env, so both lines stay in it.
doc:
	pixi run --frozen -e doc python $(CURDIR)/sh.d/make2help.py md < $(CURDIR)/Makefile > docs/MakefileTargets.md
	pixi run --frozen -e doc $(MAKE) -C docs/sphinx.d

.PHONY: cleandoc
# Clean and build documentation
cleandoc:
	pixi run --frozen -e doc $(MAKE) -C docs/sphinx.d clean
	$(MAKE) doc

.PHONY: serve-doc
# Serve built documentation at <http://localhost:8000>
serve-doc:
	pixi run --frozen -e doc python -m http.server -d docs/sphinx.d/_build/html

.PHONY: packing
# Create binary packages
packing:
	pixi run --frozen -e pack $(MAKE) -C packing

.PHONY: packing-update-containers
# Update Docker and Singularity containers used for packing
packing-update-containers:
	pixi run --frozen -e pack $(MAKE) -C packing update-containers
