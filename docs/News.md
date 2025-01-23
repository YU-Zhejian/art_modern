# News \& Release Notes

## 1.0.0

The first release of `art_modern`.

**NOTE** This version does **NOT** come with MPI support.

Changes on software function:

- [X] Supports 3 modes: `wgs`, `trans` and `templ`, similar to `pbsim3`.
- [X] Supports 3 FASTA parsers: `memory`, `htslib` and `stream`.
- [X] Supports 3 library construction methods: `se`, `pe` and `mp`.
- [X] Except FASTQ, support output in SAM and BAM format through HTSLib.
- [X] Support for masking detection was temporarily suspended.
- [X] Support for sequencers except Illumina dropped.
- [X] Support for the `aln` output format was dropped.
- [X] Built-in profiles are no longer supported. Users must specify the path to the existing profile they want to use.
- [X] Parallelization using Boost ASIO.

Changes in software implementation:

- [X] Build systems changed to CMake.
- [X] All C++ code was re-implemented in C++17 with radical removal of duplicated or unused code.
- [X] More random number generation libraries were supported.
- [X] Logging re-implemented using Boost.
- [X] Multithreading support implemented using Boost.
- [X] Largely eliminated POSIX-only routines by Boost.
- [X] Argument parser implemented in Boost.
- [X] Output writers were made asynchronous using `moodycamel::ConcurrentQueue<T>`.

## 1.0.1

Fixed miscellaneous bugs.

- Further fixed issue #2.
- More compiler versions tested; The software now supports Clang 10.0.0+, GCC 9.5.0+, and AOCC 3.2.0+.

## 1.1.0

- `--builtin_qual_file` option added back. Python 3 needed as build dependencies.
- [`BS::thread_pool`](https://github.com/bshoshany/thread-pool) added as an alternate thread pool implementation for Boost <= 1.65.
- Tested Ubuntu 18.04 x86\_64 with GCC 7.4.0, Clang 5.0.1, and Boost 1.65.1.
- Tested MacOS X Sequoia 15 with Command Line Tools for Xcode 16.2, CMake 3.31.4, and Boost 1.87.0. Fixed #3.
- Bumped bundled HTSLib to 1.21.
- Miscellaneous bug fixes.
