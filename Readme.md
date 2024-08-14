# Readme for `art_modern`

Modernized ART that is parallelized and modularized using modern C++.

**WARNING** Largely under development in internal Git hosting. The Git repository on GitHub may not reflect latest status.

## Quick Start

## Installation 

### Dependencies

- [CMake](https://cmake.org/), our building system.
- A working C and C++ compiler that supports C++ 14. Following compilers are supported:
  - [GCC](https://gcc.gnu.org/);
  - [Clang](https://clang.llvm.org/);
  - [Intel oneAPI DPCPP](https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/get-started-guide/2024-2/overview.html).
- [Boost C++ Library](https://www.boost.org/). The modules used in this project are namely:
  - **REQUIRED** [FileSystem](https://www.boost.org/doc/libs/1_85_0/libs/filesystem/);
  - **REQUIRED** [Regex](https://www.boost.org/doc/libs/1_85_0/libs/regex/);
  - **REQUIRED** [Program Options](https://www.boost.org/doc/libs/1_85_0/libs/program_options/);
  - **REQUIRED** [Thread](https://www.boost.org/doc/libs/1_85_0/libs/thread/);
  - **REQUIRED** [Log](https://www.boost.org/doc/libs/1_85_0/libs/log/);
  - **OPTIONAL** [Test](https://www.boost.org/doc/libs/1_85_0/libs/test/): For unit testing only. Can be absent for non-developers;
  - **OPTIONAL** [Timer](https://www.boost.org/doc/libs/1_85_0/libs/timer/): For displaying CPU and wall-clock time at the end of the program.
- A working [HTSLib](https://www.htslib.org/).
  - To use bundled HTSLib sources, you need to have:
    - **REQUIRED** [zlib](https://www.zlib.net/);
    - **REQUIRED** [pthread](https://www.man7.org/linux/man-pages/man7/pthreads.7.html);
    - **OPTIONAL** [libbz2](http://www.bzip.org/);
    - **OPTIONAL** [liblzma](https://tukaani.org/xz/); 
    - **OPTIONAL** [libdeflate](https://github.com/ebiggers/libdeflate);
    - See [official HTSLib documentation](https://github.com/samtools/samtools/blob/master/INSTALL) for more details.
  - To use external HTSLib, consult your system administrator.

### CMake Variables

CMake variables should be set when invoking `cmake`. For example,

```shell
cmake -DBUILD_SHARED_LIBS=ON
```

sets `BUILD_SHARED_LIBS` to `ON`.

- `BUILD_SHARED_LIBS`: Whether to build shared libraries
  - `ON` (DEFAULT): Will search for shared libraries and use dynamic linking.
  - `OFF`: Will search for static libraries and use static linking.
- `USE_HTSLIB`: Use which HTSLib implementation
  - unset (DEFAULT): Will use bundled HTSLib.
  - `hts`: Will use the HTSLib found in system.


## Usage

### Mode

### Library Construction Methods

### FASTA Parsers

## Changes Compared to Official ART Implementation

Changes on software function:

- [ ] Supports 3 modes: `wgs`, `trans` and `templ`, similar to `pbsim3`.
- [ ] Supports 3 FASTA parsers: `memory`, `htslib` and `stream`.
- [X] Supports 3 library construction methods: `se`, `pe` and `mp`.
- [X] Except FASTQ, support output in SAM and BAM format through HTSLib.
- [X] Support over masking detection was dropped.
- [X] Support over sequencers except Illumina dropped.
- [X] Support over the `aln` output format was dropped.
- [X] Built-in profiles are no longer supported. User must specify path to existing profiles.

Changes on software engineering stuff:

- Build system changed to CMake.
- All C++ code were re-implemented in C++14 with radical removal of duplicated or unused code.
- Random generator was changed from GNU Science Library (GSL) to Boost or C++ standard library.
- Logging re-implemented using Boost.
- Multithreading support implemented using Boost.
- Largely eliminated POSIX-only routines by Boost.
- Argument parser implemented in Boost.

## Acknowledgements

This simulator is based on the works of Weichun Huang <whduke@gmail.com> _et al._, under [GNU GPL v3](https://www.gnu.org/licenses/) license. The software is originally distributed [here](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) with following reference:

- W. Huang, L. Li, J. R. Myers, and G. T. Marth, _ART: a next-generation sequencing read simulator_, Bioinformatics (Oxford, England), vol. 28, no. 4, pp. 593–594, Feb. 2012, doi: [10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708).

The bundled HTSLib library used MIT License with following reference:

- J. K. Bonfield et al., _HTSlib: C library for reading/writing high-throughput sequencing data_, GigaScience, vol. 10, no. 2, p. giab007, Jan. 2021, doi: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007).

## TODO

- Support external HTSLib.
- Make it faster.
- Update the HTSLib CMake routine for setting macros like `HAVE_LIBBZ2` correct.
- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?

## FAQ

### How to split produced pair-end/mate-pair sequencing results to 2 files?
