# Installing `art_modern`

Here provides a detailed installation guide for `art_modern`.

**NOTE** All "at least" in this guide are inclusive.

## Operating System

The project assumes [x86\_64](https://en.wikipedia.org/wiki/X86-64) (aka., Intel 64, AMD64, x64, with chips manufactured by Intel, AMD, VIA, Zhaoxin, etc.) platforms. Support on other platforms (including 32-bit CPUs) are not guaranteed.

The project assumes modern GNU/Linux distributions that are still under official support. For example, Ubuntu 16.04 LTS (Xenial Xerus) had already reached its end-of-life in [Apr. 2021](https://help.ubuntu.com/community/EOL#Ubuntu_16.04_Xenial_Xerus). Other POSIX platforms like \*BSD, and patent UNIX are theoretically supported but not tested. POSIX-on-Windows platforms like [Cygwin](https://cygwin.com/), [MSYS2](https://www.msys2.org/), [MinGW](https://sourceforge.net/projects/mingw/), [MinGW-w64](https://www.mingw-w64.org/) are neither supported nor tested.

## C/C++ Compilers

### Checking C11 and C++17 Compatibility

This project requires a working C++ compiler that supports C++17 and a working C compiler that supports C11. You may test whether your compiler (GCC, for example) supports such standard using:

```shell
echo 'int main(){}' | gcc --std=c11 -x c - -o /dev/null
echo 'int main(){}' | g++ --std=c++17 -x c++ - -o /dev/null
```

**NOTE** This is a very rough test. The compiler may still fail to compile the project due to the lack of support of some specific features. In other words, if there's no error, the compiler **MIGHT** be supported. If there's an error, the compiler is definitely **NOT** supported.

**NOTE** The CMake build scripts inside this project contains a script that tests compiler compatibility which will be automatically executed when configuring the project.

For a table of the minimum compiler version that supports those versions, see also:

- [This cppreference article for C11](https://en.cppreference.com/w/c/11).
- [This cppreference article for C++17](https://en.cppreference.com/w/cpp/17).

### GCC

GNU Compiler Collections (GCC) ([homepage](https://gcc.gnu.org/)) is the most widely used compiler for GNU/Linux that provides the best compatibility and error-tolerance. At least [7.4.0](https://gcc.gnu.org/gcc-7/changes.html#GCC7.4) is required.

**NOTE** GCC supports diverse programming languages. Please ensure that your GCC installation comes with C++ support. You need at least `g++` program (Test with `g++ --version`) and a working GNU C++ Standard Library ([`libstdc++`](https://gcc.gnu.org/onlinedocs/libstdc++/)).

See also:

- GCC support over [C++17](https://gcc.gnu.org/projects/cxx-status.html#cxx17) and [C11](https://gcc.gnu.org/wiki/C11Status).
- `libstdc++` support over [C++17](https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#status.iso.2017).

### Clang

[Clang](https://clang.llvm.org/) is another popular compiler for GNU/Linux that uses Low-Level Virtual Machine ([LLVM](https://llvm.org/)) toolchain. Also, the default C++ compiler for FreeBSD and Apple Mac OS X. At least [5.0.1](https://releases.llvm.org/5.0.1/tools/clang/docs/ReleaseNotes.html) is required.

**NOTE** Clang may need GCC to work properly due to the need of the compiler runtime library (e.g., [`libgcc`/`libgcc_s`](https://gcc.gnu.org/onlinedocs/gccint/Libgcc.html) and [`libatomic`](https://gcc.gnu.org/wiki/Atomic/GCCMM)). See [this LLVM article](https://clang.llvm.org/docs/Toolchain.html) for detailed instructions on selecting GNU- or LLVM-based variants of each toolchain component for Clang.

See also:

- Clang support over [C++17](https://clang.llvm.org/cxx_status.html#cxx17) and [C11](https://clang.llvm.org/c_status.html#c11).
- `libc++` support over [C++17](https://libcxx.llvm.org/Status/Cxx17.html) if you wish to use libc++ (LLVM Standard C++ library) instead of `libstdc++` (GNU C++ Standard Library).

### Intel OneAPI DPC++/C++ Compiler

The compiler ([homepage](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html)) can build accelerated binaries on Intel CPUs. Compatibilities on the latest version of this compiler will be tested.

**NOTE** Here we refer to the LLVM-based one (with programs named `icx` and `icpx`) instead of the old Intel C++ Compiler Classic (ICC, with programs named `icc` and `icpc`).

**NOTE** CMake version higher than [3.20.0](https://cmake.org/cmake/help/latest/release/3.20.html#compilers) is required to support this compiler.

**NOTE** Distributing binaries built with this compiler may require you to comply with Intel's license and/or breaking the GPL.

### Other Compilers

Although not tested, the following compilers can also theoretically be of use:

- Intel C++ Compiler Classic (ICC) **MAY** work with legacy systems and Boost libraries.
- [NVIDIA HPC compilers](https://developer.nvidia.com/hpc-compilers) should work as all versions of this compiler support C++17.
  - **NOTE** CMake version higher than [3.20.0](https://cmake.org/cmake/help/latest/release/3.20.html#compilers) is required to support this compiler.
- [AMD Optimizing C/C++ and Fortran Compilers (AOCC)](https://www.amd.com/en/developer/aocc.html) should work as all versions of this compiler support C++17. Specifically,
  - Its first version, 3.2.0 (`AMD Clang 13.0.0 (CLANG: AOCC_3.2.0-Build#128 2021_11_12)`), was tested.
  - Its latest version, 5.0.0 (`AMD Clang 17.0.6 (CLANG: AOCC_5.0.0-Build#1377 2024_09_24)`), was tested.
- [Arm C/C++/Fortran Compiler](https://developer.arm.com/Tools%20and%20Software/Arm%20Compiler%20for%20Linux#Downloads) should work. Tested using Arm C/C++/Fortran Compiler version 24.10.1 (build number 4) (based on LLVM 19.1.0) on an ARM development board.

## Essential Tools for Building

### CMake

[CMake](https://cmake.org/) is a cross-platform build system. It is required to build the project. At least [3.17](https://cmake.org/cmake/help/latest/release/3.17.html) is required.

**NOTE** Lots of EOL distributions do not ship with a recent version of CMake. You may download CMake 3.17 in a binary form for x86\_64 GNU/Linux or Mac OS X from [CMake officially-built binaries](https://cmake.org/files/v3.17/).

### Dependencies of CMake Build System

CMake requires a [CMake Generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), which is used to perform the build. Under GNU/Linux and other POSIX systems (e.g., Mac OS X, FreeBSD), [Ninja](https://ninja-build.org/) is preferred. [GNU Make](https://www.gnu.org/software/make) is also acceptable. BSD-flavored `make` will **NOT** work.

This project requires [Python](https://www.python.org/) at least 3.7 when building the package, for embedding bundled Illumina profiles and `make help`. Python is not required at runtime.

### Linkers, Assemblers, Archivers, etc

You need either [GNU BinUtils](https://www.gnu.org/software/binutils/) or [LLVM BinUtils Replacements](https://llvm.org/docs/CommandGuide/#gnu-binutils-replacements) to perform assembling and linking. Under normal circumstances, they should come together with your compiler if you install them using your package management systems. The latter may require additional CMake variables to be set.

### C Library

This project were tested working using [GNU C Library](https://www.gnu.org/software/libc/) and [MUSL C Library](https://musl.libc.org/). Other C libraries are not tested. However, C libraries that satisfy POSIX.1-2008 and C11 should work.

**NOTE** [LLVM C Library](https://libc.llvm.org/) is neither supported nor tested.

**NOTE** Most C libraries bundles a copy of POSIX threads library ([pthread(7)](https://www.man7.org/linux/man-pages/man7/pthreads.7.html), usually named `libpthread.so` or `libpthread.a` if statically linked) and math library (Usually named `libm.so` or `libm.a`). Those libraries are required by this project.

**NOTE** Here, we assume that [`FindThread`](https://cmake.org/cmake/help/latest/module/FindThreads.html) of CMake will find pthread.

### `pkgconf` or `pkg-config`

The project may take advantage [`pkg-config`](https://www.freedesktop.org/wiki/Software/pkg-config/) or its modern replacement (**RECOMMENDED**) [`pkgconf`](https://github.com/pkgconf/pkgconf) to locate the dependencies. Installing such is **HIGHLY** recommended.

**NOTE** You may need to set up environment variable `PKG_CONFIG_PATH` to allow those tools to find the dependencies.

### Command-Line Utilities

For Apple Mac OS X, FreeBSD, Alpine Linux, and other GNU/Linux distributions with obsolete packages, please consider checking the version of the following tools:

- [GNU CoreUtils](https://www.gnu.org/software/coreutils/) (At least 8.28) for the use of `env -C` (introduced in 8.28) and `readlink -f` (introduced in 8.0).
  - **NOTE** CMake may have its own requirements on the version of GNU CoreUtils.
- [GNU Bash](https://www.gnu.org/software/bash/) (At least 4.2) for the possible use of advanced array operations.
  - **NOTE** CMake may have its own requirements on the version of GNU Bash.
- [GNU Make](https://www.gnu.org/software/make) as the BSD variant is known to be incompatible.
  - **NOTE** CMake may have its own requirements on the version of GNU Make if you use GNU Make as CMake Generator.

As your original system tools shipped with the operating system/BusyBox may **NOT** work.

## Required External Libraries

Dependencies are those libraries or tools that should be installed on your system before building the project. If you're using a personal computer with root privilege, consider installing them using your system's package manager like [APT](https://wiki.debian.org/Apt), [YUM](https://fedoraproject.org/wiki/Yum), [DNF](https://fedoraproject.org/wiki/Dnf), [`pacman`](https://wiki.archlinux.org/title/Pacman), and [Conda](https://docs.conda.io/) etc. Otherwise, contact your system administrator for where to find them or build them from source.

### Boost C++ Library

[Boost](https://www.boost.org/) is an umbrella project of diverse small modules that can be used independently. Except Boost header-only libraries, the compiled modules used in this project are:

- Essential header-only modules, including:
  - `boost/version.hpp`.
  - [Math](https://www.boost.org/doc/libs/latest/libs/math/doc/html/index.html). For calculation of the probability density function (PDF), cumulative distribution function (CDF), and inverse CDFs.
  - [Exception](https://www.boost.org/doc/libs/latest/libs/exception/doc/boost-exception.html). For better exception handling.
  - [Algorithm](https://www.boost.org/doc/libs/latest/libs/algorithm/doc/html/index.html). For simple string algorithm used in non-performance-critical situations.
- Essential modules that would be linked to the final executable:
  - [FileSystem](https://www.boost.org/doc/libs/latest/libs/filesystem/). See [this Design section](#filesystem-section) for why this module instead of C++17 `<filesystem>` is required.
  - [Program Options](https://www.boost.org/doc/libs/latest/libs/program_options/). For parsing command-line options.
  - [Log](https://www.boost.org/doc/libs/latest/libs/log/). For logging support.

**NOTE** A boost module may depend on other boost modules in either header-only or compiled form. CMake should be able to find those dependencies automatically.

Boost modules are found through `find_package(Boost ...)` command of CMake. This usually requires the presence of `BoostConfig.cmake` (Provided by Boost) or [`FindBoost.cmake`](https://cmake.org/cmake/help/latest/module/FindBoost.html) (provided by CMake) file. See [`BOOST_CONFIG_PROVIDED_BY_BOOST` CMake variable](#boost-config-provided-by-boost-section) mentioned below for details.

### zlib

[zlib](https://www.zlib.net/) performs compression and decompression bundled ART error profiles. At least [1.2.0](#v-1.2.0-section) is required. zlib is also required by [bundled HTSLib](#htslib-section).

The project will firstly try to find zlib using `pkgconf`. That usually requires the presence of `zlib.pc` file. If failed, will fall back to `libz.so`/`libz.a` with optional version suffixes.

## Optional External Libraries

The following dependencies are optional. You may choose to install them if you want to improve the performance, adaptability, or user-friendliness of the program.

(optional-boost-components)=
### Optional Boost Components

- [StackTrace](https://www.boost.org/doc/libs/latest/doc/html/stacktrace.html): For a more developer-friendly stack trace. Can be absent for non-developers.
  - See also: [configuring Boost::StackTrace](https://www.boost.org/doc/libs/latest/doc/html/stacktrace/configuration_and_build.html) for platform-specific configuration instructions for this library.
- [Test](https://www.boost.org/doc/libs/latest/libs/test/): For unit testing only. Can be absent for non-developers.
- [Timer](https://www.boost.org/doc/libs/latest/libs/timer/): For displaying CPU time, wall-clock time, and average CPU utilization at the end of the program. Can be absent if you do not care about performance.
- [Random](https://www.boost.org/doc/libs/latest/libs/random/index.html): For random number generation if you choose to use Boost random number generators (See [CMake variable `USE_RANDOM_GENERATOR`](#use-random-generator-section) below). Can be absent if you do not choose to use Boost random number generators.

If benchmarking (See [CMake variable `BUILD_ART_MODERN_BENCHMARKS`](#build-art-modern-benchmarks-section)) is required, you may also install:

- [Container](https://www.boost.org/doc/libs/latest/doc/html/container.html).
- [LockFree](https://www.boost.org/doc/libs/latest/doc/html/lockfree.html).
- [Process](https://www.boost.org/doc/libs/latest/doc/html/process.html).
- Random.

(accelerated-random-number-generators)=
### Accelerated Random Number Generators

Users on Intel/AMD CPUs are highly recommended to use [Intel OneAPI Math Kernel Library (OneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). To use this option, set [CMake variable `USE_RANDOM_GENERATOR`](#use-random-generator-section) to `ONEMKL`.

OneMKL can be found through CMake. This requires you to install the Intel OneAPI Base Toolkit, which provides `MKLConfig.cmake`. See [official docs for 2025.2](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2025-2/using-the-cmake-config-file.html) for more details.

OneMKL can also be found through `pkgconf`. For example, Debian packing of Intel MKL [`libmkl-dev`](https://packages.debian.org/sid/libmkl-dev) provides `mkl-sdl-lp64.pc`. Intel OneAPI MKL installation provides `mkl-sdl.pc`. To use MKL specified using `pkgconf` instead of CMake, see [CMake variable `FIND_RANDOM_MKL_THROUGH_PKGCONF`](#find-random-mkl-through-pkgconf-section) below.

**NOTE** OneMKL is property software. Distributing binaries built with this library may require you to comply with Intel's license and/or breaking the GPL.

[PCG](https://www.pcg-random.org/) is used as an alternative to [`std::mt19937`](https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine.html) for random number generation. File `<pcg_random.hpp>` is required. This file usually located in `/usr/include/pcg_random.hpp` if you installed Debian package [`libpcg-cpp-dev`](https://packages.debian.org/sid/all/libpcg-cpp-dev). To use this option, set [CMake variable `USE_RANDOM_GENERATOR`](#use-random-generator-section) to `SYSTEM_PCG`. Note that we also provide a bundled minified version of PCG random number generators under the name of `PCG`.

(alt-malloc-section)=
### Alternate `malloc`/`free` Implementations

Users may use either [mi-malloc](https://github.com/microsoft/mimalloc), [jemalloc](https://github.com/jemalloc/jemalloc), or [tcmalloc](https://github.com/gperftools/gperftools) to slightly improve the performance of memory allocation and deallocation and check for potential memory leaks.

- For jemalloc, `jemalloc.pc` file is required.
- For mi-malloc, `mimalloc-config.cmake` file is required.
- For tcmalloc, `libtcmalloc.pc` or `libtcmalloc_minimal.pc` is required. Under Debian GNU/Linux, this file is provided by [`libgoogle-perftools-dev`](https://packages.debian.org/sid/libgoogle-perftools-dev) package.

See also: [CMake variable `USE_MALLOC`](#use-malloc-section) below.

**NOTE** If you already compiled the binary without those libraries, you may still link them at run-time using `LD_PRELOAD`.

(mpi-section)=
### Message Passing Interface (MPI) Library

The MPI standard required in this project is MPI 1.0, which is published in 1994. Thus, most MPI libraries should work. MPI libraries from the following vendors are supported for MPI-based parallelization:

- [OpenMPI](https://www.open-mpi.org/).
- [MPICH](https://www.mpich.org/).
- [Intel MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html).

The system should theoretically support:

- [MVAPICH2](https://mvapich.cse.ohio-state.edu/).

The MPI installation with MPI C API should be locatable using CMake module [`FindMPI.cmake`](https://cmake.org/cmake/help/latest/module/FindMPI.html), which is shipped with CMake. We do **NOT** need the MPI C++ API, which is deprecated in later MPI standards.

See also: [CMake variable `WITH_MPI`](#with-mpi-section) below.

(optional-ncbi-ngs-components)=
### NCBI NGS Libraries

For enabling NCBI Short Read Archive (SRA) support in `art_profile_builder`. Requires `libncbi-ngs.so`/`libncbi-ngs.a`, with its development headers.

Available at <https://github.com/ncbi/sra-tools>.

Available since [`1.3.3`](#v-1.3.3-section).

See also: [CMake variable `WITH_NCBI_NGS`](#with-ncbi-ngs-section) below.

**NOTE** If you would like to use NCBI NGS Library, you **MUST** use bundled HTSLib due to symbol conflicts.

## Required Bundled/External Libraries

The following dependencies are bundled with the project. You do not need to install them manually. However, you may choose to use external ones if you have them installed in your system. Consult your system administrator if you do not know whether and where those libraries are installed.

**NOTE** Using bundled dependencies may introduce security vulnerabilities.

See also: [Copying.md](./Copying.md) for the licenses and versions of those bundled dependencies.

(htslib-section)=
### [HTSLib](https://www.htslib.org/)

[HTSLib](https://www.htslib.org/) is used for reading large FASTA files and generating SAM/BAM files.

See also: [CMake variable `USE_HTSLIB`](#use-htslib-section) mentioned below.

(bundled-htslib-section)=
#### Bundled

This bundled HTSLib is a trimmed down version of the original HTSLib with the following changes:

- Removed GNU Autotools \& Makefile-based building system and use CMake instead.
- Dropped support for `libcurl`.
- Removed the plugin system.
- All files unused in `htscodecs` are removed.
- Support for external HTSCodecs and OpenSSL (For calculating MD5) are removed.

To build bundled HTSLib sources, you need to have:

- **REQUIRED** zlib, which is also required by this project.
- **REQUIRED** pthread. See above section for C libraries bundling pthread.
- **HIGHLY RECOMMENDED** [libdeflate](https://github.com/ebiggers/libdeflate): This library accelerates compressed BAM output.
- **OPTIONAL** [libbz2](http://www.bzip.org/): For CRAM compression.
- **OPTIONAL** [liblzma](https://tukaani.org/xz/): For CRAM compression.

**NOTE** The CRAM format is currently not supported as an output format for `art_modern`.

See [official HTSLib documentation](https://github.com/samtools/samtools/blob/master/INSTALL) for more details.

(external-htslib-section)=
#### External

HTSLib at least [1.17](https://github.com/samtools/htslib/releases/tag/1.17) is required due to the use of `sam_flush` (1.14) and `faidx_seq_len64` (1.17).

The project will firstly try to find HTSLib using pkgconf. That usually requires the presence of `htslib.pc` file. If failed, will fall back to `lib[val].so`/`lib[val].a` with optional version suffixes where `[val]` is the value of [CMake variable `USE_HTSLIB`](#use-htslib-section).

**NOTE** Remote reference file is not supported.

### `moodycamel::ConcurrentQueue<T>`

[`moodycamel::ConcurrentQueue<T>`](https://github.com/cameron314/concurrentqueue) is used for multi-producer single-consumer output queue.

See also: [CMake variable `USE_CONCURRENT_QUEUE`](#use-concurrent-queue-section) mentioned below.

(bundled-concurrent-queue-section)=
#### Bundled

No additional dependency is required.

(external-concurrent-queue-section)=
#### External

The library is header-only, so only the path to `concurrentqueue.h` is required. At least [1.0.4](https://github.com/cameron314/concurrentqueue/releases/tag/v1.0.4) is required.

### `{fmt}`

[`{fmt}`](https://github.com/fmtlib/fmt) is used for high-speed formatting FASTA and FASTQ output.

See also: [CMake variable `USE_LIBFMT`](#use-libfmt-section) mentioned below.

(bundled-libfmt-section)=
#### Bundled

No additional dependency is required.

(external-libfmt-section)=
#### External

The project will firstly try to find `{fmt}` using pkgconf. That usually requires the presence of `fmt.pc` file. If failefd, will fall back to `lib[val].so`/`lib[val].a` with optional version suffixes, where `[val]` is the value of `USE_LIBFMT` CMake variable. At least [7.1.3](https://github.com/fmtlib/fmt/releases/tag/7.1.3) is required.

## Optional Bundled/External Dependencies

Not any.

## Required Bundled Dependencies

The required bundled dependencies of this project is `libceu`. No additional dependency is required.

## Optional Bundled Dependencies

(bs-thread-pool-section)=
### `BS::thread_pool`

[`BS::thread_pool`](https://github.com/bshoshany/thread-pool) is used as an alternative to Boost.ASIO for older Boost versions where a thread pool is not available.

No additional dependency is required.

See also: [CMake variable `USE_THREAD_PARALLEL`](#use-thread-parallel-section) mentioned below.

## CMake Variables

This project relies on diverse CMake variables that control the build behavior. If you want a specific build (e.g., with accelerated random number generation, with or without debugging information), you should set them accordingly. They should be set when invoking `cmake`. For example,

```shell
cmake -DBUILD_SHARED_LIBS=ON
```

sets `BUILD_SHARED_LIBS` to `ON`.

Following is a list of CMake variables used in this project:

### `BUILD_SHARED_LIBS`

Available since [1.0.0](#v-1.0.0-section).

This instructs CMake whether to build shared libraries. It will also affect behavior while searching for libraries.

- **`ON` (DEFAULT): Will search for shared libraries and use dynamic linking.**
- `OFF`: Will search for static libraries and use static linking.

The project should be able to be compiled into a fully static binary on [Alpine Linux](https://alpinelinux.org/) or [Void Linux](https://voidlinux.org/) with [musl libc](https://musl.libc.org/) as the standard C library. See [this blog by Li Heng](https://lh3.github.io/2014/07/12/about-static-linking) for why static linking may simplify distribution and deployment of bioinformatics software. However, this may lead to a larger binary size and security risks. See [this Debian Wiki](https://wiki.debian.org/StaticLinking) for more details.

See also: [Official CMake documentation](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html).

(cmake-build-type-section)=
### `CMAKE_BUILD_TYPE`

Available since [1.0.0](#v-1.0.0-section).

This instructs CMake to build executables/libraries with different optimization and debugging levels.

- **`Debug` (DEFAULT): For developers with debugging needs.**
  - **NOTE** Very, very slow with tons of extra checks.
- `Release`: Optimized executables/libraries without debug symbols. Used for daily use.
- `RelWithDebInfo`: Optimized executables/libraries with debug symbols. Used for profiling.

See also: [Official CMake documentation](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html).

### `Python3_EXECUTABLE`

Available since [1.1.1](#v-1.1.1-section).

Path to the Python interpreter. Required for adding bundled error profiles to the executable. Default to `python3`.

See also: [Official CMake documentation](https://cmake.org/cmake/help/latest/module/FindPython3.html).

### `CEU_CM_SHOULD_USE_NATIVE`

Available since [1.0.0](#v-1.0.0-section).

Whether to build the binaries using [`-mtune=native`](https://gcc.gnu.org/onlinedocs/gcc-14.1.0/gcc/x86-Options.html#index-march-16), if possible. This would result in faster executable with impaired portability (i.e., do not run on other machines).

- **`OFF` (DEFAULT): Will not build native executables/libraries.**
- `ON`: Will build native executables/libraries.

### `CEU_CM_SHOULD_ENABLE_TEST`

Available since [1.0.0](#v-1.0.0-section).

Whether test should be enabled.

- **Unset (DEFAULT): Set to `ON` if the CMake variable `CMAKE_BUILD_TYPE` is not `Release`, `OFF` otherwise.**
- `OFF`: Will disable test.
- `ON`: Will enable test.

(use-htslib-section)=
### `USE_HTSLIB`

Available since [1.0.0](#v-1.0.0-section).

Use which HTSLib implementation.

- **Unset (DEFAULT): Will use bundled HTSLib.** See [Bundled HTSLib](#bundled-htslib-section) for requirements.
- `hts`: Will use the HTSLib found in the system. See [External HTSLib](#external-htslib-section) for requirements.
- Any other value `[val]`: Will use the HTSLib of other names (`lib[val].so`/`lib[val].a`) found in the system. See [External HTSLib](#external-htslib-section) for requirements.

(use-random-generator-section)=
### `USE_RANDOM_GENERATOR`

Available since [1.0.0](#v-1.0.0-section).

The random number generator used.

- **`STL` (DEFAULT): Use STL random generators.**
- `PCG`: PCG random generators. This uses a minimal version of PCG random generators. Available since [1.1.1](#v-1.1.1-section).
  - Before [1.3.0](#v-1.3.0-section), this uses a bundled full-featured PCG random generators in C++.
  - After [1.3.0](#v-1.3.0-section), this uses a minimal version of PCG random generators.
- `SYSTEM_PCG`: Use PCG random generators found in the system. See See [Accelerated Random Number Generators](#accelerated-random-number-generators) for requirements. Available since [1.3.0](#v-1.3.0-section).
- `BOOST`: Use Boost random generators. See [Optional Boost Components](#optional-boost-components) for requirements.
- `ONEMKL`: Use Intel OneAPI MKL random generators. See [Accelerated Random Number Generators](#accelerated-random-number-generators) for requirements. See also [CMake variable `FIND_RANDOM_MKL_THROUGH_PKGCONF`](#find-random-mkl-through-pkgconf-section) below.

On my system (13th Gen Intel(R) Core(TM) [i7-13700H](https://www.intel.com/content/www/us/en/products/sku/232128/intel-core-i713700h-processor-24m-cache-up-to-5-00-ghz/specifications.html), Intel OneAPI BaseKit 2025.2) for generating filling 1024 random 32-bit unsigned integers 1024 times with 200 replicate, the performance is:

```text
       MKL::VSL_BRNG_SFMT19937(32 bits): gmean:        286; mean/sd:           286/3 us
         MKL::VSL_BRNG_MT19937(32 bits): gmean:        435; mean/sd:         446/132 us
               PCG::pcg32_fast(32 bits): gmean:        848; mean/sd:           848/5 us
          absl::InsecureBitGen(64 bits): gmean:      1,486; mean/sd:        1,486/33 us
        boost::random::mt19937(32 bits): gmean:      1,756; mean/sd:        1,756/19 us
                  std::mt19937(32 bits): gmean:      1,852; mean/sd:        1,853/59 us
                  GSL::mt19937(32 bits): gmean:      2,420; mean/sd:       2,422/114 us
                  absl::BitGen(64 bits): gmean:      3,543; mean/sd:       3,548/202 us
```

On another [OrangePi 3B](http://www.orangepi.org/html/hardWare/computerAndMicrocontrollers/details/Orange-Pi-3B.html) ARM development board (Rockchip RK3566 quad-core 64-bit CPU, Arm C/C++/Fortran Compiler version 24.10.1 (build number 4) (based on LLVM 19.1.0)), the performance is:

```text
               PCG::pcg32_fast(32 bits): gmean:       4,452; mean/sd:        4,452/30 us
          absl::InsecureBitGen(64 bits): gmean:      12,787; mean/sd:       12,787/66 us
        boost::random::mt19937(32 bits): gmean:      15,489; mean/sd:      15,490/117 us
                  std::mt19937(32 bits): gmean:      17,848; mean/sd:      17,848/137 us
                  GSL::mt19937(32 bits): gmean:      20,802; mean/sd:      20,802/144 us
                  absl::BitGen(64 bits): gmean:      29,104; mean/sd:      29,104/112 us
```

**NOTE** The performance may vary on different platforms.

Benchmark code available at <https://github.com/YU-Zhejian/art_modern_bench_rand>.

(use-thread-parallel-section)=
### `USE_THREAD_PARALLEL`

Available since [1.0.0](#v-1.0.0-section).

The thread-level parallelism strategy.

- **`ASIO` (DEFAULT): Will use Boost.ASIO for thread-based parallelism.**
  - **NOTE** This is only available in Boost at least 1.66.
- `BS`: Will use `BS::thread_pool`. Available since [1.1.1](#v-1.1.1-section). See [`BS::thread_pool`](#bs-thread-pool-section) for requirements.
- `NOP`: Will not use thread-based parallelism. Useful for debugging.

(boost-config-provided-by-boost-section)=
### `BOOST_CONFIG_PROVIDED_BY_BOOST`

Available since [1.1.3](#v-1.1.3-section).

Configures the behavior of CMake policy [`CMP0167`](https://cmake.org/cmake/help/latest/policy/CMP0167.html). There's usually no need to change this. You only need to set this switch to `OFF` if you have Boost earlier than 1.69 with CMake at least 3.30.

- **`ON` (DEFAULT): Will use the set the policy to `NEW`.**
- `OFF`: Will use the set the policy to `OLD`.

(use-malloc-section)=
### `USE_MALLOC`

Available since [1.1.1](#v-1.1.1-section).

Whether to use alternative high-performance `malloc`/`free` implementations like mi-malloc, jemalloc, or tcmalloc (see above). Using those implementations can improve the performance of the program but slightly increase memory consumption.

- **`AUTO` (DEFAULT): Will use jemalloc and then mi-malloc, then tcmalloc, finally minimal tcmalloc if possible.**
- `JEMALLOC`: Find and use jemalloc, and fail if not found.
- `MIMALLOC`: Find and use mi-malloc, and fail if not found.
- `TCMALLOC`: Find and use tcmalloc, and fail if not found.
- `TCMALLOC_MINIMAL`: Find and use minimal version of tcmalloc, and fail if not found.
- `NOP`: Will not use alternative `malloc`/`free` implementations. I.e., use the system-provided `malloc`/`free` implementations.

See [Alternate `malloc`/`free` Implementations](#alt-malloc-section) for requirements.

(use-libfmt-section)=
### `USE_LIBFMT`

Available since [1.1.7](#v-1.1.7-section).

Whether to use bundled `{fmt}` library for formatting strings.

- **Unset (DEFAULT): Will use bundled `{fmt}`.** See [Bundled `{fmt}`](#bundled-libfmt-section) for requirements.
- `fmt`: Will use the `{fmt}` found in the system. See [External `{fmt}`](#external-libfmt-section) for requirements.
- Any other value `[val]`: Will use the `{fmt}` of other names (`lib[val].so`) found in the system. See [External `{fmt}`](#external-libfmt-section) for requirements.

(use-concurrent-queue-section)=
### `USE_CONCURRENT_QUEUE`

Available since [1.1.7](#v-1.1.7-section).

Whether to use bundled `moodycamel::ConcurrentQueue<T>`. Specifically, where `concurrentqueue.h` is located. For example, if you use Debian GNU/Linux and installed [`libconcurrentqueue-dev`](https://packages.debian.org/sid/libconcurrentqueue-dev), you may set this variable to `/usr/include/concurrentqueue/moodycamel/` or `/usr/include/concurrentqueue/`.

- **Unset (DEFAULT): Will use bundled `moodycamel::ConcurrentQueue<T>`.** See [Bundled `moodycamel::ConcurrentQueue<T>`](#bundled-concurrent-queue-section) for requirements.
- Any value `[val]`: Will search for `moodycamel::ConcurrentQueue<T>` at including path `[val]`. See [External `moodycamel::ConcurrentQueue<T>`](#external-concurrent-queue-section) for requirements.

### `REPRODUCIBLE_BUILDS`

Available since [1.1.7](#v-1.1.7-section).

Whether to enable reproducible builds. This complies Debian policies.

- **Unset (DEFAULT): Will not enable reproducible builds.**
- `ON`: Will enable reproducible builds. All used `__DATE__` and `__TIME__` macros will be replaced with fixed values.

See also:

- [`__DATE__` and `__TIME__` in GNU C Preprocessor](https://gcc.gnu.org/onlinedocs/cpp/Standard-Predefined-Macros.html).
- [`-Wdate-time` compiler flag in Clang](https://clang.llvm.org/docs/DiagnosticsReference.html#wdate-time).
- [`-Wdate-time` compiler flag in GCC 14.3.0](https://gcc.gnu.org/onlinedocs/gcc-14.3.0/gcc/Warning-Options.html#index-Wdate-time).

(build-art-modern-benchmarks-section)=
### `BUILD_ART_MODERN_BENCHMARKS`

Available since [1.1.7](#v-1.1.7-section).

Whether to build mini benchmarks executable.

- **Unset (DEFAULT): Will not build benchmarks.**
- `ON`: Will build benchmarks.

(with-mpi-section)=
### `WITH_MPI`

Available since [1.2.0](#v-1.2.0-section).

Whether to enable MPI-based parallelization.

- **Unset (DEFAULT)/OFF: Will not enable MPI-based parallelization.**
- `ON`: Will enable MPI-based parallelization.

See [MPI Library](#mpi-section) for requirements.

(find-random-mkl-through-pkgconf-section)=
### `FIND_RANDOM_MKL_THROUGH_PKGCONF`

Available since [1.3.0](#v-1.3.0-section).

Find Intel OneAPI MKL through pkgconf. Must be specified with [CMake variable `USE_RANDOM_GENERATOR`](#use-random-generator-section) set to `ONEMKL`.

- **Unset (DEFAULT): Will not find Intel OneAPI MKL through pkgconf.**
- Any value `[val]`: Will find Intel OneAPI MKL through pkgconf with the name `[val]`. This requires the presence of `[val].pc` file.

(with-ncbi-ngs-section)=
### `WITH_NCBI_NGS`

Available since [1.3.3](#v-1.3.3-section).

Enable NCBI SRA parsing to `art_profile_builder`.

See also [NCBI NGS Libraries](#optional-ncbi-ngs-components) for requirements.

### Deprecated Options

- `USE_BTREE_MAP` was deprecated in [1.1.2](#v-1.1.2-section).
- `USE_CCACHE` was deprecated in [1.1.7](#v-1.1.7-section).
- `GSL` option of [`USE_RANDOM_GENERATOR`](#use-random-generator-section) was deprecated in [1.2.0](#v-1.2.0-section) and removed in [1.2.1](#v-1.2.1-section).
- `USE_ABSL` was deprecated and removed in [1.3.0](#v-1.3.0-section).
- `USE_QUAL_GEN` was deprecated and removed in [1.3.0](#v-1.3.0-section).
