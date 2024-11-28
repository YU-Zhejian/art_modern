# News \& Release Notes

## 1.0.0

The first release of `art_modern`.

**NOTE** This version does **NOT** come with an MPI support.

### Changes Compared to Official ART Implementation

Changes on software function:

- [X] Supports 3 modes: `wgs`, `trans` and `templ`, similar to `pbsim3`.
- [X] Supports 3 FASTA parsers: `memory`, `htslib` and `stream`.
- [X] Supports 3 library construction methods: `se`, `pe` and `mp`.
- [X] Except FASTQ, support output in SAM and BAM format through HTSLib.
- [X] Support for masking detection was dropped.
- [X] Support for sequencers except Illumina dropped.
- [X] Support for the `aln` output format was dropped.
- [X] Built-in profiles are no longer supported. Users must specify the path to the existing profile they want to use.
- [X] Parallelization using Boost ASIO.

Changes on software engineering stuff:

- Build systems changed to CMake.
- All C++ code was re-implemented in C++14 with radical removal of duplicated or unused code.
- More random number generation libraries were supported.
- Logging re-implemented using Boost.
- Multithreading support implemented using Boost.
- Largely eliminated POSIX-only routines by Boost.
- Argument parser implemented in Boost.
- FASTQ output was made asyncronous.
