# News \& Release Notes

(v-1.3.3-section)=
## 1.3.3 (2026/01/08)

- Software Engineering:
  - **EXPERIMENTAL** Added NCBI SRA support to `art_profile_builder`. Not enabled by default; Add `-DWITH_NCBI_NGS=ON` to CMake options to enable this feature.
  - Added CRAM support to `art_profile_builder`.
  - `string_dup` in bundled HTSLib renamed to `htslib_string_dup` to avoid symbol conflicts with NCBI NGS library.
  - Bumped bundled `BS::thread_pool` to [`v5.1.0`](https://github.com/bshoshany/thread-pool/releases/tag/v5.1.0).
- Implementation:
  - `art_profile_builder`: Added `--first_n_reads` to control the number of reads to be processed.
  - `art_profile_builder`: Added `--i-format` to control the input format.
  - `art_modern`: Added `--o-sam-without_tag_MD`, `--o-sam-without_tag_NM`, `--o-hl_sam-without_tag_OA`, and `--o-sam-no_qual`, to control tags of (Headless) SAM/BAM output.
- Miscellaneous bug fixes.

(v-1.3.2-section)=
## 1.3.2 (2025/12/23)

- Bug fixes:
  - The bugs affecting [`1.3.0`](#v-1.3.0-section) and [`1.3.1`](#v-1.3.1-section) should be fixed.
  - Bug involved in parsing CMake option `-DREPRODUCIBLE_BUILDS` fixed.
  - A severe bug that prevents `art_modern` from loading custom error profiles has been fixed.
- Software Engineering:
  - Bumped bundled HTSLib to [`1.23`](https://github.com/samtools/htslib/releases/tag/1.23).
- Implementation:
  - The speed of `art_profile_builder` improved by using a single-producer-multi-consumer queue. The improvement will be more obvious to SAM/BAM files and would reduce computational burden when parsing FASTQ files.
  - Queue size parameters added to `art_modern` and `art_profile_builder` to control the size of internal queues. Users may tune these parameters to achieve better performance or reduce memory usage.
  - The performance of `art_modern` under Intel OneMKL random generator further improved by generating more random numbers in bulk.
  - A seed allocation algorithm is implemented; Flag `--i-seed` added to set a fixed seed for reproduction.
- Miscellaneous bug fixes.

(v-1.3.1-section)=
## 1.3.1 (2025/12/15)

- **KNOWN BUGS**
  - Non-reproducible stucks, with no clear reason. Currently seen in both MPI and non-MPI mode. Seen in LLVM and Intel C++/DPC++ compilers. Currently, only seen in `pbsim3_transcripts`, `template` mode with unequal read1/read2 lengths specified. The author is working on fixing this issue.
- Software Engineering:
  - Update bundled `{fmt}` to [`12.1.0`](https://github.com/fmtlib/fmt/releases/tag/12.1.0).
  - A possibly unidentified bug in the CMake build system when trying to compile a static version of testing C source code while the project is being configured is fixed.
  - The minimal test code for libdeflate, zlib, libbz2, and liblzma is reimplemented in a unified style. The new behavior is more robust, as it would generate random strings, compress them, decompress them, and check whether the decompressed strings are identical to the original ones.
  - Bundled LibCEU routines can recognize more compilers \& architectures \& libraries.
- Implementation:
  - A possibly unidentified bug that would emerge when generating quality below 10 (for N bases) using Intel OneMKL random generator is fixed.
  - Several changes were made to the reporting system, increasing CPU utilization rate.
- Miscellaneous bug fixes.

(v-1.3.0-section)=
## 1.3.0 (2025/11/14)

- **BREAKING CHANGE**
  - If a 2-column (unstranded) coverage file is provided to template mode, the coverage will be interpreted as positive coverage. Prior versions would divide the coverage by 2 for each strand.
- **KNOWN BUGS**
  - The MPI mode may suffer from non-reproducible stucks, with no clear reason. This issue is not seen in non-MPI mode. The author is working on fixing this issue.
- Algorithm:
  - Introduced <https://github.com/Wunkolo/qreverse> to accelerate array reversing. Use `-DAM_NO_Q_REVERSE=ON` in CMake options to disable this feature.
  - `STL` quality generation algorithm deprecated and removed. Also deprecated CMake option `USE_QUAL_GEN`. Now `WALKER` is the only quality generation algorithm available.
  - Several low-level C++ features re-done in C to improve compilation speed.
  - `PCG` random generator: CMake option `-DUSE_RANDOM_GENERATOR=PCG` now uses a minimal PCG random generator while `-DUSE_RANDOM_GENERATOR=SYSTEM_PCG` uses an external header-only PCG random generator in C++.
- Interface:
  - **EXPERIMENTAL** `art_modern` and `art_profile_builder` now support different lengths of read 1 and 2 in paired-end simulation using `--read_len_1` and `--read_len_2` options.
    - **NOTE** For `art_modern`: Currently, different read length works for PE Template mode only. Using such on SE/MP or other modes will generate an error message when you use it.
  - `art_modern`: Some defaults added to options. `--builtin_qual_file` now default to `HiSeq2500_150bp`; `--lc` now default to `se`; `--mode` now default to `wgs`.
- Packing:
  - **EXPERIMENTAL** `-DFIND_RANDOM_MKL_THROUGH_PKGCONF` added to CMake options to find Intel OneMKL through `pkg-config`.
  - `libtcmalloc` and `libtcmalloc_minimal` supported as an alternative `malloc`/`free` implementation.
  - Google Abseil removed from dependencies. Also deprecated CMake option `USE_ABSL`.
  - `--help`, `--version` and ctest minitests removed from `debug`/`release` Makefile quick build targets and their MPI-enabled variant. Use `make testsmall`/variants instead.
- Miscellaneous bug fixes.

(v-1.2.1-section)=
## 1.2.1 (2025/11/06)

- The GNU Scientific Library (GSL) random generator was removed.
- For Intel OneMKL random generator: The bit generation routine changed to `VSL_BRNG_SFMT19937`, which is faster. Also, more random numbers are generated in bulk to reduce overhead.
- Makefile integration test targets `testbuild` and `testbuild-mpi` reimplemented in Python to make them run faster.
- `art_profile_builder` would now raise an error if the input SAM/BAM/FASTQ files are malformed.
- `art_modern`:
  - Option `--reporting_interval-job_executor` and `--reporting_interval-job_pool` added to control the reporting interval of job executor and job pool status.
  - Memory performance of `stream` FASTA parser largely improved.
- **EXPERIMENTAL** Packing: DEB package variant using OpenMPI added.
- Documentation largely revised.
- Miscellaneous bug fixes.

(v-1.2.0-section)=
## 1.2.0 (2025/10/21)

- Main repo slimmed.
  - Random number benchmark module moved to <https://github.com/YU-Zhejian/art_modern_bench_rand>.
  - Benchmark of other simulators moved to <https://github.com/YU-Zhejian/art_modern_benchmark_other_simulators>.
- The GNU Scientific Library (GSL) random generator is marked deprecated due to performance issues. They will be removed in the next release.
- **EXPERIMENTAL** Support over MPI added.
  - Currently tested MPI vendors:
    - Debian OpenMPI, `ident: 4.1.6, repo rev: v4.1.6, Sep 30, 2023`, std. 3.1.
    - Intel(R) MPI Library `2021.16 for Linux* OS`, std. 4.1.
  - MPI-related revisions:
    - Seeding of random number generators revised to avoid seed collision across different threads and processes.
    - The number of reads generated from each contig is now calculated using complicated rounding instead of flooring.
    - Makefile quick build targets `debug`, `release` have their MPI counterparts: `debug-mpi` and `release-mpi`.
    - Makefile integration test targets `testsmall`, `testsmall-release`, `testbuild` have their MPI counterparts: `testsmall-mpi`, `testsmall-release-mpi`, and `testbuild-mpi`.
    - Integration test `testbuild` revised to make it run faster.
    - CMake option `WITH_MPI` added to enable MPI support. This option is by default `OFF`.
  - **NOTE** The author currently has no access to computing clusters with MPI, so the MPI parallelization on an actual multi-node cluster may be problematic and suboptimal. Users are welcome to report bugs or tell the author how to simulate an MPI-enabled cluster using a laptop to improve the MPI support.
- Miscellaneous bug fixes.

(v-1.1.10-section)=
## 1.1.10 (2025/10/12)

- Fixed #7. In details:
  - Some build failures under Mac OS X using Apple Clang 18 have been fixed.
  - More tests added to pure-Clang/LLVM build.
- Update bundled `{fmt}` to [`12.0.0`](https://github.com/fmtlib/fmt/releases/tag/12.0.0).
- Update bundled Abseil to [`20250814.1`](https://github.com/abseil/abseil-cpp/releases/tag/20250814.1).
- Debian/Ubuntu/Alpine builder container updated to the latest versions.

(v-1.1.9-section)=
## 1.1.9 (2025/10/11)

- BAM output routines are largely accelerated by replacing string streams with pre-allocated strings.
- Implemented `art_profile_builder`, a C++ tool to build ART/`art_modern` profiles out of FASTQ or SAM/BAM files with quality scores.
- Extensively revised installation documentation.
- The simulator now supports `/dev/null` or an empty file as input.
- Static libraries removed from Alpine Linux build to reduce download size.
- HTSLib lower bound bumped to [`1.17`](https://github.com/samtools/htslib/releases/tag/1.17).
- Boost::Thread removed from required dependencies.
- Miscellaneous bug fixes.

(v-1.1.8-section)=
## 1.1.8 (2025/09/29)

- Fixed issue #5. In details:
  - On prior versions, the program will crash when trying to create simulated output in the current working directory without the `./` prefix.
  - Duplicated read IDs observed in prior versions.
  - Inconsistencies in read quality between SAM/BAM and FASTQ output were observed in prior versions.
  - Missing `/1` and `/2` suffixes in read IDs of paired-end reads observed in prior versions.
- Documentation largely revised.
- Miscellaneous bug fixes.

(v-1.1.7-section)=
## 1.1.7 (2025/09/18)

- Support for `ccache` is deprecated and removed, which also deprecated the CMake option `USE_CCACHE`.
- Update bundled Abseil to [20250814.0](https://github.com/abseil/abseil-cpp/releases/tag/20250814.0).
- Update bundled `moodycamel::ConcurrentQueue<T>` to the current latest version ([`c680721`](https://github.com/cameron314/concurrentqueue/commit/c68072129c8a5b4025122ca5a0c82ab14b30cb03)).
- Updated bundled `{fmt}` to [11.2.0](https://github.com/fmtlib/fmt/releases/tag/11.2.0).
- Some files without a clear license were removed. Unused files from bundled `{fmt}`, `moodycamel::ConcurrentQueue<T>`, and HTSLib removed.
- CMake options to use system shipped dependencies instead of bundled ones are added to comply with Debian policies. Namely, `USE_LIBFMT`, `USE_CONCURRENT_QUEUE`, `USE_ABSL`, and `REPRODUCIBLE_BUILDS`.
- Separated CMake flag that controls building of mini benchmarks to `BUILD_ART_MODERN_BENCHMARKS`.
- Miscellaneous bug fixes.

(v-1.1.6-section)=
## 1.1.6 (2025/09/12)

- A severe bug in the CMakefile of `labw_slim_htslib` was fixed.
- **EXPERIMENTAL** Debian DEB package built under Linux Mint 22 Wilma.

(v-1.1.5-section)=
## 1.1.5 (2025/09/12)

- Bumped bundled HTSLib to [`1.22.1`](https://github.com/samtools/htslib/releases/tag/1.22.1).
- Miscellaneous bug fixes.

(v-1.1.4-section)=
## 1.1.4 (2025/08/31)

- 2 environment variables, `ART_NO_LOG_DIR` and `ART_LOG_DIR`, now control the behavior of log directory creation.
- Some files without a clear license were removed.
- Add support for `cmake --install`. **NOTE** Currently, only built libraries and binaries will be installed. Documentation and header files are not included yet.
- The package is published on BioConda. See [here](https://bioconda.github.io/recipes/art_modern/README.html) for details.
- Miscellaneous bug fixes.

(v-1.1.3-section)=
## 1.1.3 (2025/03/05)

- A severe bug in built-in profiles has been fixed. Now all built-in profiles should be usable without problems. Also eliminated [-Woverlength-strings](https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html#index-Woverlength-strings) warning.
- Built-in profiles no longer being represented in base64. They're now gzip-compressed instead.
- **EXPERIMENTAL** [Sphinx](https://www.sphinx-doc.org/en/master/)-generated documentation added.
- Compiler flag `-Ofast` switched back to `-O3` for release build.
- CMake option `BOOST_CONFIG_PROVIDED_BY_BOOST` added to address newer mechanisms of Boost configuration.
- Miscellaneous bug fixes.

(v-1.1.2-section)=
## 1.1.2 (2025/02/16)

- The performance of the core simulation algorithm was improved using [Walker's Algorithm](https://doi.org/10.1145/355744.355749) on generating discrete distributions. The implementation was adapted from the C version of [GNU Scientific Library](https://www.gnu.org/software/gsl/).
- Support for B-Tree was dropped. Its performance was found to be worse than STL map in corrected benchmarks.
- The performance of `moodycamel::ConcurrentQueue<T>` was improved using producer and consumer tokens. However, since the queue is sufficiently fast without tokens, this improvement may not be significant.
- Miscellaneous bug fixes.

(v-1.1.1-section)=
## 1.1.1 (2025/02/02)

- ~~Possible build acceleration using [ccache](https://ccache.dev/) supported.~~
- Alternate `malloc`/`free` implementations like [jemalloc](https://github.com/jemalloc/jemalloc) and [mi-malloc](https://github.com/microsoft/mimalloc) are supported.
- Formatting engine of FASTQ changed to [`{fmt}`](https://github.com/fmtlib/fmt), which is slightly faster.
- FASTA output format supported.
- If the output consists only of FASTA or FASTQ, pairwise alignment will not be computed.
- The default random generator for the Intel MKL library changed from `VSL_BRNG_MT19937` to `VSL_BRNG_SFMT19937`, which is slightly faster.
- [PCG](https://www.pcg-random.org/) added as an alternative random number generator. **THIS GENERATOR MAY NOT WORK UNDER MAC OS X.**
- ~~[C++ B+ Tree](https://github.com/Kronuz/cpp-btree) added for accelerated map implementation.~~

(v-1.1.0-section)=
## 1.1.0 (2025/01/23)

- `--builtin_qual_file` option added back. Python 3 is needed as a build dependency.
- [`BS::thread_pool`](https://github.com/bshoshany/thread-pool) added as an alternate thread pool implementation for Boost less than or equal with 1.65.
- Tested Ubuntu 18.04 x86\_64 with GCC 7.4.0, Clang 5.0.1, and Boost 1.65.1.
- Tested Mac OS X Sequoia 15.2 with Command Line Tools for Xcode 16.2 (Clang 16.0.0 for target `x86_64-apple-darwin24.2.0`), CMake 3.31.4, and Boost 1.87.0. Fixed #3.
- Bumped bundled HTSLib to [`1.21`](https://github.com/samtools/htslib/releases/tag/1.21).
- Miscellaneous bug fixes.

(v-1.0.1-section)=
## 1.0.1 (2025/01/17)

Fixed miscellaneous bugs.

- Further fixed issue #2.
- More compiler versions tested; The software now supports Clang 10.0.0+, GCC 9.5.0+, and AOCC 3.2.0+.

(v-1.0.0-section)=
## 1.0.0 (2025/01/17)

The first release of `art_modern`.

**NOTE** This version does **NOT** come with MPI support.

Changes in software function:

- Supports 3 modes: `wgs`, `trans`, and `templ`, similar to `pbsim3`.
- Supports 3 FASTA parsers: `memory`, `htslib`, and `stream`.
- Supports 3 library construction methods: `se`, `pe`, and `mp`.
- Except for FASTQ, support output in SAM and BAM format through HTSLib.
- Support for masking detection was temporarily suspended.
- Support for sequencers except Illumina dropped.
- Support for the `aln` output format was dropped.
- Built-in profiles are no longer supported. Users must specify the path to the existing profile they want to use.
- Parallelization using Boost::ASIO.

Changes in software implementation:

- Build systems changed to CMake.
- All C++ code was re-implemented in C++17 with radical removal of duplicated or unused code.
- More random number generation libraries were supported.
- Logging re-implemented using Boost.
- Multithreading support implemented using Boost.
- Largely eliminated POSIX-only routines by Boost.
- Argument parser implemented in Boost.
- Output writers were made asynchronous using `moodycamel::ConcurrentQueue<T>`.
