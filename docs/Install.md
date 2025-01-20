# Installing `art_modern`

## Operating System

The project assumes x86\_64 (Intel, AMD, Zhaoxin, etc.) platforms. Other platform support is not guaranteed.

The project assumes modern GNU/Linux. That is, GNU/Linux that is under official support. For example, Ubuntu 16.04 LTS (Xenial Xerus) had already reached its end-of-life in [Apr. 2021](https://help.ubuntu.com/community/EOL#Ubuntu_16.04_Xenial_Xerus).

Other POSIX platforms like \*BSD, macOS, and patent UNIX theoretically supported but not tested. POSIX-on-Windows platforms like Cygwin, MSYS2, and MinGW are neither supported nor tested.

## Compiler Infrastructure

This project requires a working C++ compiler that supports C++17 (due to the introduction of Google ProtoBuff library) and a working C compiler that supports C11 (for bundled HTSLib).

See [here](https://en.cppreference.com/w/cpp/17) for a table of the minimum compiler version that supports C++17. You may test whether your compiler (GCC for example) supports C++17 using:

```shell
echo 'int main(){}' | g++ --std=c++17 -x c++ - -o /dev/null
```

If there's no error, the compiler is supported. You may also test your compiler using [`MericLuc/Cpp17-Features-tests`](https://github.com/MericLuc/Cpp17-Features-tests/) or other C++17 test suites.

### [GCC](https://gcc.gnu.org/)

The most widely used compiler for GNU/Linux that provides the best compatibility and error-tolerance.

- **NOTE** GCC supports diverse programming languages. Please ensure that your GCC installation comes with C++ support. You need at least `g++` program (Test with `g++ --version`) and a working GNU C++ Standard Library (libstdc++).
- The minimal version of GCC tested is GCC 7.4.0.
- Reference:
  - GCC support over [C++17](https://gcc.gnu.org/projects/cxx-status.html#cxx17) and [C11](https://gcc.gnu.org/wiki/C11Status).
  - libstdc++ support over [C++17](https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#status.iso.2017).

### [Clang](https://clang.llvm.org/)

Another popular compiler for GNU/Linux that uses [LLVM](https://llvm.org/) toolchain. Also, the default C++ compiler for FreeBSD and Apple macOS.

- **NOTE** Clang may need GCC to work properly due to the need of the compiler runtime library ([`libgcc`/`libgcc_s`](https://gcc.gnu.org/onlinedocs/gccint/Libgcc.html) and [`libatomic`](https://gcc.gnu.org/wiki/Atomic/GCCMM)). See [here](https://clang.llvm.org/docs/Toolchain.html) for detailed instructions on selecting GNU- or LLVM-based variants of each toolchain component for Clang.
- The minimal version of Clang tested is Clang 5.0.1.
- However, earlier Clang versions may be supported with earlier operating systems, GCC, and Boost library.
- Reference:
  - Clang support over [C++17](https://clang.llvm.org/cxx_status.html#cxx17) and [C11](https://clang.llvm.org/c_status.html#c11).
  - Libc++ support over [C++17](https://libcxx.llvm.org/Status/Cxx17.html) if you wish to use libc++ (LLVM Standard C++ library) instead of libstdc++ (GNU C++ Standard Library).
  - [LLVM C Library](https://libc.llvm.org/) is neither supported nor tested.

### [Intel oneAPI DPC++/C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html)

For accelerated binaries on Intel CPUs. Note that here we refer to the LLVM-based one (with programs named `icx` and `icpx`) instead of the old Intel C++ Compiler Classic (abbr., ICC, with programs named `icc` and `icpc`).

### Other Compilers

Although not tested, the following compilers can also theoretically be of use:

- Intel C++ Compiler Classic (ICC) **MAY** work with legacy systems and Boost libraries.
- [NVidia HPC compilers](https://developer.nvidia.com/hpc-compilers).
  - All versions of this compiler support C++17.
- [AMD Optimizing C/C++ and Fortran Compilers (AOCC)](https://www.amd.com/en/developer/aocc.html).
  - All (>=3.2.0) versions of this compiler support C++17. Specifcically,
    - Its first version, 3.2.0 (`AMD Clang 13.0.0 (CLANG: AOCC_3.2.0-Build#128 2021_11_12)`), was tested.
    - Its latest version, 5.0.0 (`AMD Clang 17.0.6 (CLANG: AOCC_5.0.0-Build#1377 2024_09_24)`), was tested.

### Essential Tools for Building

You need either [GNU BinUtils](https://www.gnu.org/software/binutils/) or [LLVM BinUtils Replacements](https://llvm.org/docs/CommandGuide/#gnu-binutils-replacements) to perform assembling and linking. The latter may require additional CMake variables to be set.

For C library, this project works on [GNU C Library](https://www.gnu.org/software/libc/) and [MUSL C Library](https://musl.libc.org/). Other C libraries are not tested.

## Using CMake Building System

[CMake](https://cmake.org/) 3.17 or above is recommended to build this project, although this project **MAY** work with CMake >= 3.10. That further requires a [CMake Generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), which is used to perform the build. Under GNU/Linux and other POSIX systems (e.g., macOS, FreeBSD), using [Ninja](https://ninja-build.org/) is preferred. [GNU Make](https://www.gnu.org/software/make) is also acceptable.

The CMake modules used in this project further requires [Python](https://www.python.org/) >= 3.7 and a POSIX-compiliant shell (e.g., [Dash](http://gondor.apana.org.au/~herbert/dash/), [Bash](https://www.gnu.org/software/bash/), etc.) for several text processing functions.

This project relies on diverse CMake variables that control the build behavior. If you want a specific build (e.g., with accelerated random number generation, with or without debugging information), you should set them accordingly. They should be set when invoking `cmake`. For example,

```shell
cmake -DBUILD_SHARED_LIBS=ON
```

sets `BUILD_SHARED_LIBS` to `ON`.

## Dependencies

Dependencies are those libraries or tools that should be installed on your system before building the project. If you're using a personal computer with root privilege, consider installing them using your system's package manager like [APT](https://wiki.debian.org/Apt), YUM/DNF, PacMan, etc. Otherwise, contact your system administrator for where to find them or build them from source.

### [Boost C++ Library](https://www.boost.org/)

This is an umbrella project of diverse small modules that can be used independently. The modules used in this project are namely:

- **REQUIRED** [FileSystem](https://www.boost.org/doc/libs/1_85_0/libs/filesystem/).
- **REQUIRED** [Regex](https://www.boost.org/doc/libs/1_85_0/libs/regex/).
- **REQUIRED** [Program Options](https://www.boost.org/doc/libs/1_85_0/libs/program_options/).
- **REQUIRED** [Thread](https://www.boost.org/doc/libs/1_85_0/libs/thread/).
- **REQUIRED** [Log](https://www.boost.org/doc/libs/1_85_0/libs/log/).
- **REQUIRED** [StackTrace](https://www.boost.org/doc/libs/1_85_0/doc/html/stacktrace.html).
- **OPTIONAL** [Test](https://www.boost.org/doc/libs/1_85_0/libs/test/): For unit testing only. Can be absent for non-developers.
- **OPTIONAL** [Timer](https://www.boost.org/doc/libs/1_85_0/libs/timer/): For displaying CPU and wall-clock time at the end of the program. Can be absent if you do not care about performance.
- **OPTIONAL** `stacktrace_backtrace`: For a more developer-friendly stack trace. See [here](https://www.boost.org/doc/libs/1_85_0/doc/html/stacktrace/configuration_and_build.html) for details. Can be absent for non-developers.

- **NOTE** Clang <=9 with Boost >= 1.76 may raise bug when using headers from `boost/math` as boost may misidentify Clang as GCC that does not support C++11. See [here](https://www.boost.org/doc/libs/1_87_0/libs/math/doc/html/math_toolkit/history2.html#math_toolkit.history2.math_3_0_0_boost_1_76) for more details.
- **NOTE** Boost <=1.65 does not contain `boost/asio/thread_pool.hpp` or `boost/asio/post.hpp`. See [here](https://www.boost.org/doc/libs/1_66_0/doc/html/boost_asio/reference/thread_pool.html) for its first introduction.

### [HTSLib](https://www.htslib.org/)

You may either use the one bundled with the project or an external one that had already been installed inside your system. To use bundled HTSLib sources, you need to have:

- **REQUIRED** [zlib](https://www.zlib.net/).
- **REQUIRED** [pthread](https://www.man7.org/linux/man-pages/man7/pthreads.7.html).
- **OPTIONAL** [libbz2](http://www.bzip.org/): For CRAM compression, which is now not supported.
- **OPTIONAL** [liblzma](https://tukaani.org/xz/): For CRAM compression, which is now not supported.
- **OPTIONAL** [libdeflate](https://github.com/ebiggers/libdeflate): This library accelerates compressed BAM output.

See [official HTSLib documentation](https://github.com/samtools/samtools/blob/master/INSTALL) for more details. See also `USE_HTSLIB` CMake variable mentioned below.

To use external HTSLib, consult your system administrator. Those libraries usually named `libhts.so` with optional version suffixes.

### Other Optional Libraries

- Various libraries for accelerating random number generation. Including:
  - **OPTIONAL** [Intel OneAPI Math Kernel Library (OneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html).
  - **OPTIONAL** [GNU Science Library (GSL)](https://www.gnu.org/software/gsl/).

## Reference on CMake Variables

### [`BUILD_SHARED_LIBS`](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html)

This instructs CMake whether to build shared libraries. It will also affect behavior while searching for libraries.

- **`ON` (DEFAULT): Will search for shared libraries and use dynamic linking.**
- `OFF`: Will search for static libraries and use static linking.

### [`CMAKE_BUILD_TYPE`](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html)

This instructs CMake to build executables/libraries with different optimization and debugging levels.

- **`Debug` (DEFAULT): For developers with debugging needs.**
- `Release`: Optimized executables/libraries without debug symbols. Used for daily use.
- `RelWithDebInfo`: Optimized executables/libraries with debug symbols. Used for profiling.

### `Python3_EXECUTABLE`

Path to [Python](https://www.python.org/) >= 3.7. Default to `python3`.

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
- `BOOST`: Use Boost random generators.
- `GSL`: Use GSL random generators. This is used in the original ART.
- `ONEMKL`: Use Intel OneAPI MKL random generators. Highly recommended on Intel CPUs.

### `USE_THREAD_PARALLEL`

The thread-level parallelism strategy.

- **`ASIO` (DEFAULT): Will use Boost.ASIO for thread-based parallelism.**
  - **NOTE** This is only available in Boost >= 1.66.
- `BS`: Will use [`BS::thread_pool`](https://github.com/bshoshany/thread-pool).
- `NOP`: Will not use thread-based parallelism. Useful for debugging.

### `BOOST_CONFIG_PROVIDED_BY_BOOST`

Configures the behaviour of CMake policy [`CMP0167`](https://cmake.org/cmake/help/latest/policy/CMP0167.html).

- **`ON` (DEFAULT): Will use the set the policy to `NEW`.**
- `OFF`: Will use the set the policy to `OLD`.

There's ususally no need to change this. You only need to set this switch to `OFF` if you have Boost < 1.70 with CMake >= 3.30.

## Appendix

```shell
CLANG_TARNAME="clang+llvm-10.0.0-x86_64-linux-gnu-ubuntu-18.04"
CLANG_DIR="${HOME}/opt/${CLANG_TARNAME}/"

env -C "${HOME}/opt/" wget https://github.com/llvm/llvm-project/releases/download/llvmorg-10.0.0/${CLANG_TARNAME}.tar.xz
env -C "${HOME}/opt/" tar xvJf "${CLANG_TARNAME}.tar.xz"

make clean release \
    PATH="${CLANG_DIR}/bin/:${PATH}" \
    CMAKE_FLAGS="-DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
```

Test script for ICC 2023.2.4, its latest release:

```shell
. /opt/intel/oneapi/compiler/2023.2.4/env/vars.sh
make clean release \
    PATH="${CLANG_DIR}/bin/:${PATH}" \
    CMAKE_FLAGS="-DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc"
```
