# Developing \& Contributing

## Common Development-Oriented Tasks

### Create Development Environment

Except [CMake](https://cmake.org) and other dependencies specified in [Install](Install.md), the development scripts also require other dependencies. Install [Conda](https://docs.conda.io/en/latest/) (Or [Mamba](https://mamba.readthedocs.io/en/latest/)/[MicroMamba](https://mamba.readthedocs.io/en/latest/micromamba.html)), and then execute:

```shell
conda env create -f env/art_modern.yml
```

### Testing

Some individual modules can be tested using CMake CTest system. Execute `ctest` in your building directory to test them.

Run `make testsmall` or `make testsmall-release` to run the integration tests using executables produced in `make debug` and `make release`.

### Profiling

We have `./profile.sh ${PROFILER}` to perform profiling of the software with profiler `${PROFILER}`. Profilers may include:

- `intel-advisor` for Intel Advisor. Intel compilers with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.
- `intel-vtune` for Intel VTune Profiler. Intel compilers with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.
- `nsys` for NVIDIA Nsight Systems. NVIDIA compiler with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.
- `valgrind` for Valgrind Callgrind. GCC with `Debug` mode will be used since optimization may produce instructions that are not supported by Valgrind.
- `amd-uprof` for AMD uProf. AMD compilers with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.

### Building Documentations

Create the development environment and make the documentation through:

```shell
conda activate art_modern
make doc
```

And the built documentations (HTML and PDF) should be in `doc/sphinx.d/_build/html` and `doc/sphinx.d/_build/latex` respectively. Note that for PDF output, you may need [latexmk](https://www.ctan.org/pkg/latexmk) and a working up-to-date [LaTeX](https://www.latex-project.org) distribution (e.g., [TeXLive](https://www.tug.org/texlive/), [MiKTeX](https://miktex.org/), [MacTeX](https://tug.org/mactex/)) installed.

## Release Cycle

- Feature freeze: All new features merged into the `develop` branch.
- Code freeze: All **production** code changes (including bug fixes) merged into the `develop` branch. Run `make fmt` to format the code. Ensure that integration test passes.
- Integration test: `make testbuild` and `make testbuild-mpi` are executed using GCC, pure LLVM/Clang, and Intel OneAPI C++/DPC++ toolchain to ensure compatibility. Ensure `make packing` passes without errors.
- Documentation freeze: All documentation changes merged into the `develop` branch. Note:
  - Update `News.md`.
  - Update version number in `CMakeLists.txt`.
  - Ensure `make cleandoc` passes without errors.
- Release: The `develop` branch is merged into the main branch, tagged with a new version number.
- Generate artifact: Run `make packing` and `make cleandoc`.

## Get Engaged

### Issues

You're welcome to submit issues on GitHub if you've encountered any problems while using this software. You're recommended to take advantage of the bug report templates, which will help us to understand your problem better. You may also submit issues if you're unclear about the docs or have some ideas for improving the software itself.

### Pull Requests (PRs)

You're welcomed to send pull requests (PRs) to this project using the standard fork-and-pull-request workflow. However, please send an issue first to discuss the changes you're going to make (except spelling \& punctuation \& grammar issues). Otherwise, we may not be able to accept your PR.

Before you send a PR, please make sure that:

- You've run `make fmt` to format the code.
- All CTest passes.
- `make clean testsmall-release testsmall-release-mpi` passes. This test usually requires 6 to 8 minutes.
- `make clean testsmall testsmall-mpi` passes. This test usually requires 12 to 18 minutes. Here, the additional tests that are not enabled in release mode will be activated.
- `make testbuild testbuild-mpi` passes. This test usually requires 4 to 6 hours.

You may also:

- You used Valgrind to check for memory leaks.

## Scripts for Testing

### Pure LLVM/Clang Toolchain

A LLVM toolchain file is provided in `sh.d/toolchain/host-llvm/llvm-toolchain.cmake` for users who want to use Clang/LLVM toolchain. This toolchain uses LLVM `libc++` as C++ standard library and LLVM `ld.lld` as linker.

**NOTE** The Boost library shipped through your system may be compiled with GNU C++ ABI, which is not compatible with LLVM `libc++`. You may need to build Boost from source using Clang/LLVM toolchain. All C libraries do not have this issue since C ABI is generally compatible across different compilers under GNU/Linux.

```shell
LD_LIBRARY_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${LD_LIBRARY_PATH:-}" \
    LD_RUN_PATH="${HOME}/opt/boost-1.89.0-clang/lib/:${LD_RUN_PATH:-}"\
    PKG_CONFIG_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/:${PKG_CONFIG_PATH:-}" \
    CMAKE_TOOLCHAIN_FILE="$(pwd)/sh.d/toolchain/host-llvm/llvm-toolchain.cmake" \
    make testbuild testbuild-mpi \
    CMAKE_FLAGS='-DBoost_DIR=${HOME}/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/'
```

Given that `{fmt}` and Boost are installed in `${HOME}/opt/fmt-12.0.0-clang` and `${HOME}/opt/boost-1.89.0-clang` respectively.

### Intel DPCPP Toolchain

```shell
. /opt/intel/oneapi/setvars.sh
make testbuild testbuild-mpi \
    CMAKE_FLAGS='-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx'
```

## Miscellaneous Developer-Oriented Documentation

```{toctree}
MakefileTargets.md
TODO.md
CODE_OF_CONDUCT.md
```
