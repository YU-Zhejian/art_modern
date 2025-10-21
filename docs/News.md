# News \& Release Notes

## 1.2.0 (Ongoing)

- Random number benchmark module moved to <https://github.com/YU-Zhejian/art_modern_bench_rand>.
- The GNU Science Library (GSL) random generator is marked deprecated due to performance issues. They will be removed in the next release.
- **EXPERIMENTAL** Support over MPI added.
  - Currently tested MPI vendors:
    - Debian OpenMPI, ident: 4.1.6, repo rev: v4.1.6, Sep 30, 2023.
  - MPI-related revisions:
    - Seeding of random number generators revised to advoid seed colission accross different threads and processes.
    - The number of reads generated from each contig is now calculated using complicated rounding instead of flooring, which may increase the total number of reads by a small amount but overall makes the capturing more precise.
      - In details, the number of reads generated at positive and negative strands will be adjusted by 1 (SE) or 2 (PE/MP) to make the number of generated bases and the number of required bases as close as possible.
      - For example, consider generating 125-nt reads 5.0 positive and 5.0 negative depth for a contig of length 225.
        - In SE mode, 9 reads will be generated on positive and negative strands.
        - In PE/MP mode, 10 and 8 reads will be generated on positive and negative strands respectively.
    - Makefile quick build targets `debug`, `release` have their MPI counterparts: `debug-mpi` and `release-mpi`.
    - Makefile integration test targets `testsmall`, `testsmall-release`, `testbuild` have their MPI counterparts: `testsmall-mpi`, `testsmall-release-mpi`, and `testbuild-mpi`.
    - Integration test `testbuild` revised to make it run faster.
    - CMake option `WITH_MPI` added to enable MPI support. This option is by default `OFF`.
  - **NOTE** The author currently have no access to computing clusters with MPI, so the MPI parallelization on an actual multi-node cluster may be problematic and suboptimal. Users are welcomed to report bugs or tell the author how to simulate MPI-enabled cluster using a laptop to improve the MPI support.
- Miscellaneous bug fixes.

## 1.1.10 (2025/10/12)

- Fixed #7. In details:
  - Some build failure under Mac OS X using Apple Clang 18 fixed.
  - More tests added to pure-Clang/LLVM build.
- Update bundled `{fmt}` to [`12.0.0`](https://github.com/fmtlib/fmt/releases/tag/12.0.0).
- Update bundled Abseil to [`20250814.1`](https://github.com/abseil/abseil-cpp/releases/tag/20250814.1).
- Debian/Ubuntu/Alpine builder container updated to the latest versions.

## 1.1.9 (2025/10/11)

- BAM output routines are largely accelerated by replacing string streams with pre-allocated strings.
- Implemented `art_profile_builder`, a C++ tool to build ART/`art_modern` profiles out of FASTQ or SAM/BAM files with quality scores.
- Extensively revised documentation for installation.
- The simulator now supports `/dev/null` or an empty file as input.
- Static libraries removed from Alpine Linux build to reduce download size.
- HTSLib lower bound bumped to 1.17.
- Boost.Thread removed from required dependencies.
- Miscellaneous bug fixes.

## 1.1.8 (2025/09/29)

- Fixed issue #5. In details:
  - On prior versions, the program will crash when trying to create simulated output in the current working directory without `./` prefix.
  - Duplicated read IDs observed in prior versions.
  - Inconsistencies in read quality between SAM/BAM and FASTQ output observed in prior versions.
  - Missing `/1` and `/2` suffixes in read IDs of paired-end reads observed in prior versions.
- Documentation largely revised.
- Miscellaneous bug fixes.

## 1.1.7 (2025/09/18)

- Support over `ccache` deprecated, which also deprecated CMake option `USE_CCACHE`.
- Update bundled Abseil to [20250814.0](https://github.com/abseil/abseil-cpp/releases/tag/20250814.0).
- Update bundled `moodycamel::ConcurrentQueue<T>` to the current latest version ([`c680721`](https://github.com/cameron314/concurrentqueue/commit/c68072129c8a5b4025122ca5a0c82ab14b30cb03)).
- Updated bundled `{fmt}` to [11.2.0](https://github.com/fmtlib/fmt/releases/tag/11.2.0).
- Some files without clear license were removed. Unused files from bundled `{fmt}`, `moodycamel::ConcurrentQueue<T>`, and HTSLib removed.
- CMake options to use system shipped dependencies instead of bundled ones are added to comply with Debian policies. Namely, `USE_LIBFMT`, `USE_CONCURRENT_QUEUE`, `USE_ABSL`, and `REPRODUCIBLE_BUILDS`.
- Separated CMake flag that controls building of mini benchmarks to `BUILD_ART_MODERN_BENCHMARKS`.
- Miscellaneous bug fixes.

## 1.1.6 (2025/09/12)

- A severe bug in CMakefiles of `labw_slim_htslib` fixed.
- **EXPERIMENTAL** Debian DEB package built under Linux Mint 22 Wilma.

## 1.1.5 (2025/09/12)

- Bumped bundled HTSLib to 1.22.1.
- Miscellaneous bug fixes.

## 1.1.4 (2025/08/31)

- 2 environment variables, `ART_NO_LOG_DIR` and `ART_LOG_DIR` now controllers the behavior of log directory creation.
- Some files without clear license were removed.
- Add support for `cmake --install`. **NOTE** Currently, only built libraries and binaries will be installed. Documentation and header files are not included yet.
- The package is published at BioConda. See [here](https://bioconda.github.io/recipes/art_modern/README.html) for details.
- Miscellaneous bug fixes.

## 1.1.3 (2025/03/05)

- A severe bug in builtin profiles fixed. Now all builtin profiles should be usable without problems. Also eliminated [-Woverlength-strings](https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html#index-Woverlength-strings) warning.
- Builtin profiles no longer being represented in base64. They're now gzip-compressed instead.
- **EXPERIMENTAL** [Sphinx](https://www.sphinx-doc.org/en/master/)-generated documentation added.
- Compiler flag `-Ofast` switched back to `-O3` for release build.
- CMake option `BOOST_CONFIG_PROVIDED_BY_BOOST` added to address newer mechanisms of Boost configuration.
- Miscellaneous bug fixes.

## 1.1.2 (2025/02/16)

- The performance of the core simulation algorithm was improved using [Walker's Algorithm](https://doi.org/10.1145/355744.355749) on generating discrete distributions. The implementation was adapted from the C version of [GNU Science Library](https://www.gnu.org/software/gsl/).
- Support over B-Tree was dropped. Its performance was found worse than STL map in corrected benchmarks.
- The performance of `moodycamel::ConcurrentQueue<T>` was improved using producer and consumer tokens. However, since the queue is sufficiently fast without tokens, this improvement may not be significant.
- Miscellaneous bug fixes.

## 1.1.1 (2025/02/02)

- ~~Possible build acceleration using [ccache](https://ccache.dev/) supported.~~
- Alternate `malloc`/`free` implementations like [jemalloc](https://github.com/jemalloc/jemalloc) and [mi-malloc](https://github.com/microsoft/mimalloc) supported.
- Formatting engine of FASTQ changed to [`{fmt}`](https://github.com/fmtlib/fmt), which is slightly faster.
- FASTA output format supported.
- If the output consists only FASTA or FASTQ, pairwise alignment will not be computed.
- The default random generator for the Intel MKL library changed from `VSL_BRNG_MT19937` to `VSL_BRNG_SFMT19937`, which is slightly faster.
- [PCG](https://www.pcg-random.org/) added as an alternative random number generator. **THIS GENERATOR MAY NOT WORK UNDER MAC OS X.**
- ~~[C++ B+ Tree](https://github.com/Kronuz/cpp-btree) added for accelerated map implementation.~~

## 1.1.0 (2025/01/23)

- `--builtin_qual_file` option added back. Python 3 needed as build dependencies.
- [`BS::thread_pool`](https://github.com/bshoshany/thread-pool) added as an alternate thread pool implementation for Boost <= 1.65.
- Tested Ubuntu 18.04 x86\_64 with GCC 7.4.0, Clang 5.0.1, and Boost 1.65.1.
- Tested Mac OS X Sequoia 15.2 with Command Line Tools for Xcode 16.2 (Clang 16.0.0 for target `x86_64-apple-darwin24.2.0`), CMake 3.31.4, and Boost 1.87.0. Fixed #3.
- Bumped bundled HTSLib to 1.21.
- Miscellaneous bug fixes.

## 1.0.1 (2025/01/17)

Fixed miscellaneous bugs.

- Further fixed issue #2.
- More compiler versions tested; The software now supports Clang 10.0.0+, GCC 9.5.0+, and AOCC 3.2.0+.

## 1.0.0 (2025/01/17)

The first release of `art_modern`.

**NOTE** This version does **NOT** come with MPI support.

Changes on software function:

- Supports 3 modes: `wgs`, `trans` and `templ`, similar to `pbsim3`.
- Supports 3 FASTA parsers: `memory`, `htslib` and `stream`.
- Supports 3 library construction methods: `se`, `pe` and `mp`.
- Except FASTQ, support output in SAM and BAM format through HTSLib.
- Support for masking detection was temporarily suspended.
- Support for sequencers except Illumina dropped.
- Support for the `aln` output format was dropped.
- Built-in profiles are no longer supported. Users must specify the path to the existing profile they want to use.
- Parallelization using Boost ASIO.

Changes in software implementation:

- Build systems changed to CMake.
- All C++ code was re-implemented in C++17 with radical removal of duplicated or unused code.
- More random number generation libraries were supported.
- Logging re-implemented using Boost.
- Multithreading support implemented using Boost.
- Largely eliminated POSIX-only routines by Boost.
- Argument parser implemented in Boost.
- Output writers were made asynchronous using `moodycamel::ConcurrentQueue<T>`.
