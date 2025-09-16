# Installing `art_modern`

## Operating System

The project assumes x86\_64 (Intel, AMD, Zhaoxin, etc.) platforms. Other platform support is not guaranteed.

The project assumes modern GNU/Linux distributions that are still under official support. For example, Ubuntu 16.04 LTS (Xenial Xerus) had already reached its end-of-life in [Apr. 2021](https://help.ubuntu.com/community/EOL#Ubuntu_16.04_Xenial_Xerus).

Other POSIX platforms like \*BSD, and patent UNIX are theoretically supported but not tested. POSIX-on-Windows platforms like [Cygwin](https://cygwin.com/), [MSYS2](https://www.msys2.org/), [MinGW](https://sourceforge.net/projects/mingw/), [MinGW-w64](https://www.mingw-w64.org/) are neither supported nor tested.

This project used bundled source code of [Google Abseil](https://abseil.io/), whose requirements are available [here](https://github.com/google/oss-policies-info/blob/main/foundational-cxx-support-matrix.md).

## Compiler Infrastructure

This project requires a working C++ compiler that supports C++17 and a working C compiler that supports C11 (for bundled HTSLib). That also includes the C++ standard library, compiler runtime library, and miscellaneous tools like linker and assembler.

### Checking C11 and C++17 Compatibility

See [here](https://en.cppreference.com/w/c/11) and [here](https://en.cppreference.com/w/cpp/17) for a table of the minimum compiler version that supports C11 and C++17. You may test whether your compiler (GCC, for example) using:

```shell
echo 'int main(){}' | gcc --std=c11 -x c - -o /dev/null
echo 'int main(){}' | g++ --std=c++17 -x c++ - -o /dev/null
```

If there's no error, the compiler is supported. The CMake build scripts inside this project also contains a script that tests compiler compatibility using codes from [GNU AutoConf Archive](https://www.gnu.org/software/autoconf-archive/).

### [GCC](https://gcc.gnu.org/)

The most widely used compiler for GNU/Linux that provides the best compatibility and error-tolerance. The minimal version of GCC tested is GCC 7.4.0.

**NOTE** GCC supports diverse programming languages. Please ensure that your GCC installation comes with C++ support. You need at least `g++` program (Test with `g++ --version`) and a working GNU C++ Standard Library ([`libstdc++`](https://gcc.gnu.org/onlinedocs/libstdc++/)).

Reference:

- GCC support over [C++17](https://gcc.gnu.org/projects/cxx-status.html#cxx17) and [C11](https://gcc.gnu.org/wiki/C11Status).
- `libstdc++` support over [C++17](https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#status.iso.2017).

### [Clang](https://clang.llvm.org/)

Another popular compiler for GNU/Linux that uses [LLVM](https://llvm.org/) toolchain. Also, the default C++ compiler for FreeBSD and Apple Mac OS X. The minimal version of Clang tested is Clang 5.0.1.

**NOTE** Clang may need GCC to work properly due to the need of the compiler runtime library ([`libgcc`/`libgcc_s`](https://gcc.gnu.org/onlinedocs/gccint/Libgcc.html) and [`libatomic`](https://gcc.gnu.org/wiki/Atomic/GCCMM)). See [here](https://clang.llvm.org/docs/Toolchain.html) for detailed instructions on selecting GNU- or LLVM-based variants of each toolchain component for Clang.

Reference:

- Clang support over [C++17](https://clang.llvm.org/cxx_status.html#cxx17) and [C11](https://clang.llvm.org/c_status.html#c11).
- `libc++ `support over [C++17](https://libcxx.llvm.org/Status/Cxx17.html) if you wish to use libc++ (LLVM Standard C++ library) instead of `libstdc++` (GNU C++ Standard Library).

**NOTE** [LLVM C Library](https://libc.llvm.org/) is neither supported nor tested.

### [Intel oneAPI DPC++/C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html)

For accelerated binaries on Intel CPUs. Note that here we refer to the LLVM-based one (with programs named `icx` and `icpx`) instead of the old Intel C++ Compiler Classic (abbr., ICC, with programs named `icc` and `icpc`). Only the latest version is tested.

**NOTE** Distributing binaries built with this compiler may require you to comply with Intel's license.

### Other Compilers

Although not tested, the following compilers can also theoretically be of use:

- Intel C++ Compiler Classic (ICC) **MAY** work with legacy systems and Boost libraries.
- [NVidia HPC compilers](https://developer.nvidia.com/hpc-compilers) should work as all versions of this compiler support C++17. However, old CMake may not be able to recognize the compiler.
- [AMD Optimizing C/C++ and Fortran Compilers (AOCC)](https://www.amd.com/en/developer/aocc.html).
  - All (>=3.2.0) versions of this compiler support C++17. Specifically,
    - Its first version, 3.2.0 (`AMD Clang 13.0.0 (CLANG: AOCC_3.2.0-Build#128 2021_11_12)`), was tested.
    - Its latest version, 5.0.0 (`AMD Clang 17.0.6 (CLANG: AOCC_5.0.0-Build#1377 2024_09_24)`), was tested.

## Essential Tools for Building

### [CMake](https://cmake.org/)

CMake 3.17 or above is recommended to build this project. That further requires a [CMake Generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), which is used to perform the build. Under GNU/Linux and other POSIX systems (e.g., Mac OS X, FreeBSD), using [Ninja](https://ninja-build.org/) is preferred. [GNU Make](https://www.gnu.org/software/make) is also acceptable.

The CMake modules used in this project further require [Python](https://www.python.org/) >= 3.7 and a POSIX-compliant shell (e.g., [Dash](http://gondor.apana.org.au/~herbert/dash/), [Bash](https://www.gnu.org/software/bash/), etc.) for several text processing functions.

### Linkers, Assemblers, Archivers, etc

You need either [GNU BinUtils](https://www.gnu.org/software/binutils/) or [LLVM BinUtils Replacements](https://llvm.org/docs/CommandGuide/#gnu-binutils-replacements) to perform assembling and linking. Under normal circumstances, they should come together with your compiler if you install them using your package management systems. The latter may require additional CMake variables to be set.

### C Library

This project works on [GNU C Library](https://www.gnu.org/software/libc/) and [MUSL C Library](https://musl.libc.org/). Other C libraries are not tested. However, C libraries that satisfy POSIX.1-2008 should work.

### [`pkgconf`](https://github.com/pkgconf/pkgconf)  or [`pkg-config`](https://www.freedesktop.org/wiki/Software/pkg-config/)

The project may take advantage of pkgconfor pkg-config to locate the dependencies. Installing such is **HIGHLY** recommended.

### Command-Line Ultilities

For Apple Mac OS X, FreeBSD, Alpine Linux, and other GNU/Linux distributions with obsolete packages, please consider check the version of the following tools:

- [GNU CoreUtils](https://www.gnu.org/software/coreutils/) (At least 8.28) for the use of `env -C` (introduced in 8.28) and `readlink -f` (introduced in 8.0).
  - **NOTE** CMake may have its own requirements on the version of GNU CoreUtils.
- [GNU Bash](https://www.gnu.org/software/bash/) (At least 4.2) for the possible use of advanced array operations.
  - **NOTE** CMake may have its own requirements on the version of GNU Bash.
- [GNU Make](https://www.gnu.org/software/make) (Lowest version TBD) as the BSD variant is known to be incompatible.
  - **NOTE** CMake may have its own requirements on the version of GNU Make if you use GNU Make as CMake Generator.

as your original system tools shipped with the operating system/BusyBox may **NOT** work.

## Library Dependencies

Dependencies are those libraries or tools that should be installed on your system before building the project. If you're using a personal computer with root privilege, consider installing them using your system's package manager like [APT](https://wiki.debian.org/Apt), [YUM](https://fedoraproject.org/wiki/Yum), [Dnf](https://fedoraproject.org/wiki/Dnf), [`pacman`](https://wiki.archlinux.org/title/Pacman), etc. Otherwise, contact your system administrator for where to find them or build them from source.

### [Boost C++ Library](https://www.boost.org/)

This is an umbrella project of diverse small modules that can be used independently. Except Boost header-only libraries, the compiled modules used in this project are namely:

- **REQUIRED** [FileSystem](https://www.boost.org/doc/libs/1_85_0/libs/filesystem/).
- **REQUIRED** [Program Options](https://www.boost.org/doc/libs/1_85_0/libs/program_options/).
- **REQUIRED** [Thread](https://www.boost.org/doc/libs/1_85_0/libs/thread/).
- **REQUIRED** [Log](https://www.boost.org/doc/libs/1_85_0/libs/log/).
- **OPTIONAL** [StackTrace](https://www.boost.org/doc/libs/1_85_0/doc/html/stacktrace.html): For a more developer-friendly stack trace. See [here](https://www.boost.org/doc/libs/1_85_0/doc/html/stacktrace/configuration_and_build.html) for details. Can be absent for non-developers.
- **OPTIONAL** [Test](https://www.boost.org/doc/libs/1_85_0/libs/test/): For unit testing only. Can be absent for non-developers.
- **OPTIONAL** [Timer](https://www.boost.org/doc/libs/1_85_0/libs/timer/): For displaying CPU and wall-clock time at the end of the program. Can be absent if you do not care about performance.

### [HTSLib](https://www.htslib.org/)

You may either use the one bundled with the project or an external one that had already been installed inside your system. To use bundled HTSLib sources, you need to have:

- **REQUIRED** [zlib](https://www.zlib.net/).
- **REQUIRED** [pthread](https://www.man7.org/linux/man-pages/man7/pthreads.7.html).
- **OPTIONAL** [libbz2](http://www.bzip.org/): For CRAM compression. Note that CRAM format is currently not supported as an output format for `art_modern`.
- **OPTIONAL** [liblzma](https://tukaani.org/xz/): For CRAM compression, which is now not supported.
- **OPTIONAL** [libdeflate](https://github.com/ebiggers/libdeflate): This library accelerates compressed BAM output. **HIGHLY RECOMMENDED**.

See [official HTSLib documentation](https://github.com/samtools/samtools/blob/master/INSTALL) for more details. See also `USE_HTSLIB` CMake variable mentioned below.

To use external HTSLib, consult your system administrator. Those libraries usually named `libhts.so`/`libhts.a` with optional version suffixes.

### [zlib](https://www.zlib.net/)

For compression and decompression bundled ART error profiles. At least 1.2.0 (launched 2003/03/09) is required.

### Other Optional Libraries

- Various libraries for accelerating random number generation. Including:
  - **OPTIONAL** [Intel OneAPI Math Kernel Library (OneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). Highly recommended if you are using Intel CPUs.
    - **NOTE** OneMKL is available for Intel CPUs only.
    - **NOTE** This library is propertary software.
  - **OPTIONAL** [GNU Science Library (GSL)](https://www.gnu.org/software/gsl/).
- Alternate `malloc`/`free` implementations that may improve performance. Including:
  - **OPTIONAL** [mi-malloc](https://github.com/microsoft/mimalloc).
  - **OPTIONAL** [jemalloc](https://github.com/jemalloc/jemalloc).

## CMake Variables

This project relies on diverse CMake variables that control the build behavior. If you want a specific build (e.g., with accelerated random number generation, with or without debugging information), you should set them accordingly. They should be set when invoking `cmake`. For example,

```shell
cmake -DBUILD_SHARED_LIBS=ON
```

sets `BUILD_SHARED_LIBS` to `ON`.

Following is a list of CMake variables used in this project:

### [`BUILD_SHARED_LIBS`](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html)

This instructs CMake whether to build shared libraries. It will also affect behavior while searching for libraries.

- **`ON` (DEFAULT): Will search for shared libraries and use dynamic linking.**
- `OFF`: Will search for static libraries and use static linking.

The project should be able to be compiled into a fully static binary on [Alpine Linux](https://alpinelinux.org/) or [Void Linux](https://voidlinux.org/) with [musl libc](https://musl.libc.org/) as the standard C library. See [this blog by Li Heng](https://lh3.github.io/2014/07/12/about-static-linking) for why static linking may simplify distribution and deployment of bioinformatics software. However, this may lead to a larger binary size and security risks. See [this Debian Wiki](https://wiki.debian.org/StaticLinking) for more details.

### [`CMAKE_BUILD_TYPE`](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html)

This instructs CMake to build executables/libraries with different optimization and debugging levels.

- **`Debug` (DEFAULT): For developers with debugging needs.**
  - **NOTE** Very, very slow with tons of extra checks.
- `Release`: Optimized executables/libraries without debug symbols. Used for daily use.
- `RelWithDebInfo`: Optimized executables/libraries with debug symbols. Used for profiling.

### [`Python3_EXECUTABLE`](https://cmake.org/cmake/help/latest/module/FindPython3.html)

Path to [Python](https://www.python.org/) >= 3.7. Default to `python3`.

Available since 1.1.1.

### `CEU_CM_SHOULD_USE_NATIVE`

Whether to build the binaries using [`-mtune=native`](https://gcc.gnu.org/onlinedocs/gcc-14.1.0/gcc/x86-Options.html#index-march-16), if possible. This would result in faster executable with impaired portability (i.e., do not run on other machines).

- **`OFF` (DEFAULT): Will not build native executables/libraries.**
- `ON`: Will build native executables/libraries.

### `CEU_CM_SHOULD_ENABLE_TEST`

Whether test should be enabled.

- **unset (DEFAULT): Set to `ON` if the CMake variable `CMAKE_BUILD_TYPE` is not `Release`, `OFF` otherwise.**
- `OFF`: Will disable test.
- `ON`: Will enable test.

### `USE_HTSLIB`

Use which HTSLib implementation.

- **unset (DEFAULT): Will use bundled HTSLib.**
- `hts`: Will use the HTSLib (`libhts.so`) found in the system.
- Any other value `[val]`: Will use the HTSLib of other names (`lib[val].so`) found in the system.

### `USE_RANDOM_GENERATOR`

The random number generator used.

- **`STL` (DEFAULT): Use STL random generators.**
- `PCG`: [PCG](https://www.pcg-random.org/) random generators. Available since 1.1.1.
  - **NOTE** Experimental.
  - **NOTE** This generator would fail on Mac OS X due to the lack of `cxxabi.h`.
- `BOOST`: Use Boost random generators.
- `GSL`: Use GSL random generators.
  - This is used in the original ART.
  - **NOTE** Slow on Intel CPUs.
- `ONEMKL`: Use Intel OneAPI MKL random generators.
  - **NOTE** Highly recommended on Intel CPUs.

On my system for generating 1M random bits for 20K times:

```text
PCG::pcg32(0, 4294967295): 2,840 us
MKL::VSL_BRNG_SFMT19937 (32 bits): 6,472 us
MKL::VSL_BRNG_MT19937 (32 bits): 8,672 us
boost::random::mt19937(0, 4294967295): 12,669 us
std::mt19937(0, 4294967295): 19,837 us
absl::InsecureBitGen(0, 18446744073709551615): 30,928 us
GSL::mt19937(0, 4294967295): 48,060 us
```

### `USE_THREAD_PARALLEL`

The thread-level parallelism strategy.

- **`ASIO` (DEFAULT): Will use Boost.ASIO for thread-based parallelism.**
  - **NOTE** This is only available in Boost >= 1.66.
- `BS`: Will use [`BS::thread_pool`](https://github.com/bshoshany/thread-pool). Available since 1.1.1.
- `NOP`: Will not use thread-based parallelism. Useful for debugging.

### `BOOST_CONFIG_PROVIDED_BY_BOOST`

Configures the behavior of CMake policy [`CMP0167`](https://cmake.org/cmake/help/latest/policy/CMP0167.html).

- **`ON` (DEFAULT): Will use the set the policy to `NEW`.**
- `OFF`: Will use the set the policy to `OLD`.

There's usually no need to change this. You only need to set this switch to `OFF` if you have Boost < 1.70 with CMake >= 3.30.

### `USE_QUAL_GEN`

Available since 1.1.2.

The quality generation algorithm.

- **`WALKER` (DEFAULT): Use [Walker's Algorithm](https://doi.org/10.1145/355744.355749) to accelerate quality synthesis.**
- `STL`: Use binary search algorithm implemented in `std::map`. This is identical to the original ART.

### `USE_MALLOC`

Available since 1.1.1.

Whether to use alternative high-performance `malloc`/`free` implementations like mi-malloc or jemalloc. Using those implementations can improve the performance of the program but slightly increase memory consumption.

- **`AUTO` (DEFAULT): Will use jemalloc and then mi-malloc if possible.**
- `JEMALLOC`: Find and use jemalloc, and fail if not found.
- `MIMALLOC`: Find and use mi-malloc, and fail if not found.
- `NOP`: Will not use alternative `malloc`/`free` implementations. I.e., use the system-provided `malloc`/`free` implementations.

### `USE_LIBFMT`

Available since 1.1.7.

Whether to use bundled `{fmt}` library for formatting strings.

- **unset (DEFAULT): Will use bundled `{fmt}`.**
- `fmt`  : Will use the `{fmt}` (`libfmt.so`) found in the system.
- Any other value `[val]`: Will use the `{fmt}` of other names (`lib[val].so`) found in the system.

### `USE_CONCURRENT_QUEUE`

Available since 1.1.7.

Whether to use bundled `moodycamel::ConcurrentQueue`.

- **unset (DEFAULT): Will use bundled `moodycamel::ConcurrentQueue`.**
- Any value `[val]`: Will search for MoodyCamel ConcurrentQueue at including path `[val]`. For example, if you use Debian GNU/Linux and intalled [`libconcurrentqueue-dev`](https://packages.debian.org/sid/libconcurrentqueue-dev), you may set this variable to `/usr/include/concurrentqueue/moodycamel/`.

### `USE_ABSL`

Available since 1.1.7.

Whether to use bundled Abseil library.

- **unset (DEFAULT): Will use bundled Abseil.**
- Any value `[val]`: Will use system Abseil found by the CMake module `abslConfig.cmake`, which is shipped with official Abseil libraries.

### Deprecated Options

- `USE_BTREE_MAP` was deprecated in 1.1.2.
- `USE_CCACHE` was deprecated in 1.1.7.

## Platform-Specific Building Instructions

### Debian, Ubuntu, Linux Mint, or Other Debian-Based Distributions

Install the following packages using APT.

```shell
apt install -y build-essential g++ binutils libboost-all-dev libz-dev make apt-file python3 cmake
```

And then you may build this project using CMake.

**NOTE** If the CMake version is too old, you may need to install a newer version from the official repository.

### Alpine Linux

For building a static binary, you may need to install the following packages using APK:

```shell
apk add g++ binutils boost-dev zlib-dev make python3 cmake zlib-static icu-static coreutils boost1.84-static
```

**NOTE** Please install the correct version of Boost static library.

**NOTE** `coreutils` is **MANDATORY** -- Those shipped with BusyBox will **NOT** work.

### Apple Mac OS X

Install Xcode command-line tools using:

```shell
xcode-select --install
```

Alternatively, you may get the latest version [here](https://developer.apple.com/download/all). An Apple account is required.

Download CMake from [here](https://cmake.org/download/). You should modify PATH in `~/.bashrc` or `~/.zshrc` to include the `bin` directory of the CMake installation. If you install the DMG image, it will commonly be located in `/Applications/CMake.app/Contents/bin/cmake`.

Download Boost source code from [here](https://www.boost.org/users/download/). Extract it. Build it from source using:

```shell
# Assume you got Boost 1.87.0.
# Assume you've extracted the source code in your Downloads version
cd Downloads/boost_1_87_0
./bootstrap.sh
./b2
```

And then you may use CMake to build this project through:

```shell
mkdir -p build
cd build
# Set -DBoost_DIR accordingly.
# Older CMake may have different behaviour.
cmake .. -DBoost_DIR=/Users/USERNAME/Downloads/boost_1_87_0/stage/lib/cmake/Boost-1.87.0
```

You may also set up dependencies using [Conda](https://docs.conda.io), [MacPorts](https://www.macports.org/) or [HomeBrew](https://brew.sh).

## Building Documentations

Create the development environment and make the documentation through:

```shell
conda activate art_modern
make doc
```

And the built documentations (HTML and PDF) should be in `doc/sphinx.d/_build/html` and `doc/sphinx.d/_build/latex` respectively. Note that for PDF output, you may need [latexmk](https://www.ctan.org/pkg/latexmk) and a working up-to-date [LaTeX](https://www.latex-project.org) distribution (e.g., [TeXLive](https://www.tug.org/texlive/), [MiKTeX](https://miktex.org/), [MacTeX](https://tug.org/mactex/)) installed.

This step can be executed under Microsoft Windows if you have some GNU Make installed. If you use MikTeX, it may complain on not finding a `perl`. But you can have it installed through [Conda](https://docs.conda.io) anyway through:

```shell
conda install -c conda-forge perl
```

## Known Incompatibilities

- Clang < 10 with Boost > 1.77 may raise bug when using headers from `boost/math` as boost may misidentify Clang as GCC that does not support C++11. See [here](https://www.boost.org/doc/libs/1_87_0/libs/math/doc/html/math_toolkit/history2.html#math_toolkit.history2.math_3_0_0_boost_1_76) for more details.
- Boost < 1.66 does not contain `boost/asio/thread_pool.hpp` or `boost/asio/post.hpp`. See [here](https://www.boost.org/doc/libs/1_66_0/doc/html/boost_asio/reference/thread_pool.html) for its first introduction. Under this circumstance, try `-DUSE_THREAD_PARALLEL=BS` in CMake options (introduced below) as an alternative.
- GCC 13.3.0 on Haiku OS hrev58590 may generate a kernel panic that jam the entire system.
- GCC would fail on Debian GNU/Hurd.
