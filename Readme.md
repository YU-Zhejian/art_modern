# Readme for `art_modern`

Modernized ART that is parallelized and modularized using modern C++.

Changes:

- Supports 3 modes: `wgs`, `trans` and `templ`, similar to `pbsim3`.
- Supports 3 library construction methods: `se`, `pe` and `mp`.
- Supports Illumina only. Support over the `aln` output format was dropped.
- Build system changed to CMake.
- All C++ code were re-implemented in C++14 with radical removal of duplicated or unused code.
- Random generator was changed from GNU Science Library (GSL) to Boost or C++ standard library.
- Logging re-implemented using Boost.
- Multithreading support implemented using Boost.
- Largely eliminated POSIX-only routines by Boost.
- Argument parser implemented in Boost.
- Built-in profiles are no longer supported. User must specify path to existing profiles.

This simulator is based on the works of Weichun Huang <whduke@gmail.com>, under [GNU GPL v3](https://www.gnu.org/licenses/) license.
