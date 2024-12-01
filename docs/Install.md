# Installing `art_modern`

## Compiler

This project requires a working C++ compiler that supports C++17 (due to the introduction of Google ProtoBuff library) and a working C compiler that supports C11 (for bundled HTSLib).

The following compilers are tested and supported:

- [GCC](https://gcc.gnu.org/): The most widely-used compiler for GNU/Linux that provides the best compatibility and error-tolerance.
  - See [here](https://gcc.gnu.org/projects/cxx-status.html#cxx17) for support over C++17 in GCC.
- [Clang](https://clang.llvm.org/): Another popular compiler for GNU/Linux. Also, the default C++ compiler for FreeBSD and Apple macOS.
  - See [here](https://clang.llvm.org/cxx_status.html#cxx17) for support over C++17 in Clang.
- [Intel oneAPI DPC++/C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html): For accelerated binaries on Intel CPUs.

Although not tested, the following compilers can also theoretically be of use:

- Intel C++ Compiler Classic (ICC).
- [NVidia HPC compilers](https://developer.nvidia.com/hpc-compilers).
- [AMD Optimizing C/C++ and Fortran Compilers (AOCC)](https://www.amd.com/en/developer/aocc.html).

See [here](https://en.cppreference.com/w/cpp/17) for a table of the minimum compiler version that supports C++17. You may test whether your compiler (GCC for example) supports C++17 using:

```shell
echo 'int main(){}' | g++ --std=c++17 -x c++ - -o /dev/null
```

If there's no error, the compiler is supported.

## Using CMake Building System

[CMake](https://cmake.org/) 3.17 or above is required to build this project. That further requires a [CMake Generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), which is used to perform the build. Under GNU/Linux and other POSIX systems (e.g., macOS, FreeBSD), using [Ninja](https://ninja-build.org/) is preferred. [GNU Make](https://www.gnu.org/software/make) is also acceptable.

This project relies on diverse CMake variables that control the build behaviour. If you want a specific build (e.g., with accelerated random number generation, with or without debugging information), you should set them accordingly. They should be set when invoking `cmake`. For example,

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

This instructs CMake whether to build shared libraries. It will also affect behaviour while searching for libraries.

- **`ON` (DEFAULT): Will search for shared libraries and use dynamic linking.**
- `OFF`: Will search for static libraries and use static linking.

### [`CMAKE_BUILD_TYPE`](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html)

This instructs CMake to build executables/libraries with different optimisation and debugging levels.

- **`Debug` (DEFAULT): For developers with debugging needs.**
- `Release`: Optimized executables/libraries without debug symbols. Used for daily use.
- `RelWithDebInfo`: Optimized executables/libraries with debug symbols. Used for profiling.

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

- **`ASIO` (DEFAULT): Will use Boost ASIO for thread-based parallelism.**
- `NOP`: Will not use thread-based parallelism. Useful for debugging.
