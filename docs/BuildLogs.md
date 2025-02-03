# Build Logs

Here contains some build logs for unsupported platforms. If you insist on building `art_modern` on those platforms, they may be of help.

## Ubuntu 18.04 chroot

Git commit: `8ccdf5f57ce81993cecebe116007c3b73bb28c3c`.

Host: Linux Mint 12 with `6.8.0-51-generic` x86\_64 kernel.

External CMake (cmake-3.31.4-linux-x86_64 from official website) was used.

Installed packages:

| Package            | Version                 |
|--------------------|-------------------------|
| `build-essential`  | 12.4ubuntu1             |
| `cmake`            | 3.31.4                  |
| `g++`              | 7.4.0-1ubuntu2.3        |
| `clang++-5`        | 1:5.0.1-4               |
| `binutils`         | 2.30-21ubuntu1~18.04.9  |
| `libboost-all-dev` | 1.65.1.0ubuntu1         |
| `make`             | 4.1-9.1ubuntu1          |
| `libstdc++6`       | 8.4.0-1ubuntu1~18.04    |
| `libc6`            | 2.27-3ubuntu1.6         |
| `python3`          | 3.6.7-1~18.04           |

`art_modern --version` output:

```text
ART: 2.5.8, ART_MODERN: 1.1.0
ART_MODERN_LINK_LIBS: Boost::filesystem;Boost::regex;Boost::program_options;Boost::thread;Boost::log_setup;Boost::log;labw_slim_htslib;libceu
USING HTSLib: labw_slim_htslib, ver. 1.21
	Features: build=Makefile libcurl=no S3=no GCS=no libdeflate=no lzma=no bzip2=no plugins=no htscodecs=1.6.1
	CC: /usr/bin/cc
	CPPFLAGS: 
	CFLAGS: -O0;-pg;-g3;-Og;-pedantic;-Wextra;-Wall
	LDFLAGS: 
	HTSlib URL scheme handlers present:
		built-in: file, data, preload
		crypt4gh-needed: crypt4gh
		mem: mem
BOOST: 1.65.1
GSL: not used
MKL: not used
MPI: not used
Protobuf: not used
OpenMP: not used
SIMDE: N/A
	with ISE: MMX SSE2
Compile-time C std.: ver. unknown ISO C (__STDC_VERSION__ undefined)
Compile-time C++ std.: ver. 17 (201703)
Compiled at: Jan 19 2025, 14:40:12
Compiler Identification:
	GCC compatible version number: 7.5.0
	__VERSION__ string: 7.5.0
Compile-time C Types max, min, etc. limits:
	char           (1 size):                       -128 ->                  +127
	schar          (1 size):                       -128 ->                  +127
	uchar          (1 size):                         +0 ->                  +255
	size_t         (8 size):                         +0 ->  18446744073709551615
	ptrdiff_t      (8 size):       -9223372036854775808 ->  +9223372036854775807
	short          (2 size): [...]
	ushort         (2 size): [...]
	int            (4 size): [...]
	uint           (4 size): [...]
	long           (8 size): [...]
	ulong          (8 size): [...]
	llong          (8 size): [...]
	ullong         (8 size): [...]
	bool           (1 size): [...]
Compile-time OS info:
	PRIMIARY OS='GNU/Linux'
	POSIX.1 Version: 200809
	POSIX.2 Version: 200809
	Single UNIX Specification (SUS) Version: 700
Run-time OS info:
	POSIX UTSINFO:
		sysname=Linux
		nodename=yuzj-HP-ENVY-x360-Convertible-15-dr1xxx
		release=6.8.0-51-generic
		version=#52-Ubuntu SMP PREEMPT_DYNAMIC Thu Dec  5 13:09:44 UTC 2024
		machine=x86_64
```

```text
[...]
Compiler Identification:
	Clang compatible version number: 5.0.1 (5.0.1 (tags/RELEASE_501/final))
	GCC compatible version number: 4.2.1
	__VERSION__ string: 4.2.1 Compatible Clang 5.0.1 (tags/RELEASE_501/final)
[...]
```

Other notes:

- Threading disabled. Boost 1.65.1 does not provide `boost::asio::thread_pool`, so [`BS::thread_pool`](https://github.com/bshoshany/thread-pool) is used instead.
- Final chroot size: 2.0G.

## FreeBSD 14.2-RELEASE VM

Git commit: `1703e45833d9b4481b9d45b155407d5ec23eb3c9`.

FreeBSD: 14.2-RELEASE with x86\_64 kernel.

Installed packages:

| Package            | Version                 |
|--------------------|-------------------------|
| `cmake`            | 3.31.4                  |
| `binutils`         | 2.43.1.1  |
| `cmake`            | 3.31.3 |
| `boost-all`| 1.86.0         |
| `gmake`             | 4.4.1          |
| `python311`          | 3.11.11           |

`art_modern --version` output:

```text
ART: 2.5.8, ART_MODERN: 1.1.2
On git commit: (1703e45833d9b4481b9d45b155407d5ec23eb3c9) Mon, 3 Feb 2025 20:00:28 +0800
ART_MODERN_LINK_LIBS: Boost::filesystem;Boost::program_options;Boost::thread;Boost::log_setup;Boost::log;Boost::timer;Boost::stacktrace_basic;labw_slim_htslib;slim_libceu;slim_libfmt;Threads::Threads
USING HTSLib: labw_slim_htslib, ver. 1.21
	Features: build=Makefile libcurl=no S3=no GCS=no libdeflate=yes lzma=yes bzip2=yes plugins=no htscodecs=1.6.1
	CC: /usr/bin/cc
	CPPFLAGS: 
	CFLAGS: -DHAVE_LIBLZMA;-DHAVE_LIBDEFLATE;-DHAVE_LIBBZ2;-DHAVE_DRAND48;-Ofast;-g0;-w;-mtune=native;-march=native
	LDFLAGS: CEU_CM_EFL::libm_shared;CEU_CM_EFL::libz_shared;CEU_CM_EFL::libbz2_shared;CEU_CM_EFL::libdeflate_shared;CEU_CM_EFL::liblzma_shared
	HTSlib URL scheme handlers present:
		built-in: file, data, preload
		crypt4gh-needed: crypt4gh
		mem: mem
{fmt}: 7.1.3
BOOST: 1.86.0
GSL: not used
MKL: not used
MPI: not used
Protobuf: not used
OpenMP: not used
SIMDE: N/A
	with ISE: MMX SSE2 AVX2
BS::thread_pool: 5.0.0
*malloc: not used
Compile-time C std.: ver. unknown ISO C (__STDC_VERSION__ undefined)
Compile-time C++ std.: ver. 17 (201703)
Compiled at: Feb  3 2025, 20:01:06
Compiler Identification:
	Clang compatible version number: 18.1.6 (18.1.6 (https://github.com/llvm/llvm-project.git llvmorg-18.1.6-0-g1118c2e05e67))
	GCC compatible version number: 4.2.1
	__VERSION__ string: FreeBSD Clang 18.1.6 (https://github.com/llvm/llvm-project.git llvmorg-18.1.6-0-g1118c2e05e67)
Compile-time C Types max, min, etc. limits:
	char           (1 size):                       -128 ->                  +127
	schar          (1 size):                       -128 ->                  +127
	uchar          (1 size):                         +0 ->                  +255
	size_t         (8 size):                         +0 ->  18446744073709551615
	ptrdiff_t      (8 size):       -9223372036854775808 ->  +9223372036854775807
	short          (2 size): [...]
	ushort         (2 size): [...]
	int            (4 size): [...]
	uint           (4 size): [...]
	long           (8 size): [...]
	ulong          (8 size): [...]
	llong          (8 size): [...]
	ullong         (8 size): [...]
	bool           (1 size): [...]
Compile-time OS info:
	PRIMIARY OS='Free BSD'
	POSIX.1 Version: 200112
	POSIX.2 Version: 199212
	Single UNIX Specification (SUS) Version: unknown
Run-time OS info:
	POSIX UTSINFO:
		sysname=FreeBSD
		nodename=ylbsd
		release=14.2-RELEASE
		version=FreeBSD 14.2-RELEASE releng/14.2-n269506-c8918d6c7412 GENERIC
		machine=amd64
```
