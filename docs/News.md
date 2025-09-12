# News \& Release Notes

## 1.1.6 (2025/09/12)

- A severe bug in CMakefiles of labw_slim_htslib fixed.

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
- Miscellaneous bug fixes.

Bundled files:

- `build_rel_with_dbg_alpine-x86_64.zip`: Static linked libraries and executable binaries built under x86\_64 Alpine Linux. Should work on most x86\_64 Linux distributions.
- `art_modern.pdf`: Documentation in PDF format.
- `art_modern.html.zip`: Documentation in HTML format.

## 1.1.2 (2025/02/16)

- The performance of the core simulation algorithm was improved using [Walker's Algorithm](https://doi.org/10.1145/355744.355749) on generating discrete distributions. The implementation was adapted from the C version of [GNU Science Library](https://www.gnu.org/software/gsl/).
- Support over B-Tree was dropped. Its performance was found worse than STL map in corrected benchmarks.
- The performance of MoodyCamel queue was improved using producer and consumer tokens. However, since the queue is sufficiently fast without tokens, this improvement may not be significant.
- Miscellaneous bug fixes.

Bundled files:

- `build_rel_with_dbg_alpine-x86_64.zip`: Static linked libraries and executable binaries built under x86\_64 Alpine Linux. Should work on most x86\_64 Linux distributions.

## 1.1.1 (2025/02/02)

- Possible build acceleration using [ccache](https://ccache.dev/) supported.
- Alternate `malloc`/`free` implementations like [jemalloc](https://github.com/jemalloc/jemalloc) and [mi-malloc](https://github.com/microsoft/mimalloc) supported.
- Formatting engine of FASTQ changed to [`{fmt}`](https://github.com/fmtlib/fmt), which is slightly faster.
- FASTA output format supported.
- If the output consists only FASTA or FASTQ, pairwise alignment will not be computed.
- The default random generator for the Intel MKL library changed from `VSL_BRNG_MT19937` to `VSL_BRNG_SFMT19937`, which is slightly faster.
- [PCG](https://www.pcg-random.org/) added as an alternative random number generator. **THIS GENERATOR MAY NOT WORK UNDER MAC OS X.**
- ~~[C++ B+ Tree](https://github.com/Kronuz/cpp-btree) added for accelerated map implementation.~~

Bundled files:

- `art_modern_alpine`: Static linked binary built under x86\_64 Alpine Linux. Should work on most x86\_64 Linux distributions.

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

