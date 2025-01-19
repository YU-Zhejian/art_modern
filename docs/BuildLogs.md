# Build Logs

## Ubuntu 18.04 chroot

Installed packages:

| Package | Version |
|---------|---------|
| `build-essential` | 12.4ubuntu1 |
| `cmake` | 3.26.4 |
| `g++` | 7.4.0-1ubuntu2.3 |
| `binutils` | 2.30-21ubuntu1~18.04.9 |
| `cmake` | 3.10.2-1ubuntu2.18.04.2 |
| `libboost-all-dev` | 1.65.1.0ubuntu1 |
| `make` | 4.1-9.1ubuntu1 |
| `libstdc++6` | 8.4.0-1ubuntu1~18.04 |

`art_modern --version` output:

```text
[2025-01-19 15:10:58.504369] : YuZJ Modified ART_Illumina (art_modern v. 1.1.0)
[2025-01-19 15:10:58.504639] : Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)
[2025-01-19 15:10:58.504665] : Originally written by: Weichun Huang <whduke@gmail.com>
[2025-01-19 15:10:58.504684] : Modified by: YU Zhejian <Zhejianyu@intl.zju.edu.cn>
[2025-01-19 15:10:58.504703] : Debugging functions enabled.
[2025-01-19 15:10:58.504860] : MPI not found! Cross-node parallelism disabled.
[2025-01-19 15:10:58.505499] : ARGS: /root/src/build/art_modern --version
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
[2025-01-19 15:10:58.506351] : EXIT
```
