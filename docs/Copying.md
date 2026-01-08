# Copyright

This project uses GPL v3 License, a copy of which is available at [License.md](../License.md).

---

This project is based on the source code of [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) by [Weichun Huang](mailto:whduke@gmail.com) _et al._, under the [GNU GPL v3](https://www.gnu.org/licenses/) license.

---

This project uses code from the following projects (in no particular order). All code used is either in the public domain or under a license that is [compatible with GPL v3](https://www.gnu.org/licenses/license-list.en.html#GPLCompatibleLicenses).

## HTSlib by Genome Research Ltd

Available from <https://github.com/samtools/htslib>.

Affected files:

- `/deps/labw_slim_htslib/cram/**`, under the MIT/Expat License.
- `/deps/labw_slim_htslib/**`, under The Modified-BSD License.

**NOTE** This project uses the source code retrieved from [GitHub 1.23 Release Tarball](https://github.com/samtools/htslib/releases/download/1.23/htslib-1.23.tar.bz2).

## `gitignore` by GitHub

Available from <https://github.com/github/gitignore>.

Affected files:

- `/.gitignore`

With CC0 1.0 Universal License.

**NOTE** Added in [1.0.0](#v-1.0.0-section).

## `libceu` by YU Zhejian

Affected files:

- `/deps/cmake_collections/**`
- `/deps/slim_libceu/**`
  
With MIT License.

**NOTE** Added in [1.0.0](#v-1.0.0-section).

**NOTE** This project was dead. New projects should consider [Boost::Predef](https://www.boost.org/doc/libs/1_87_0/libs/predef/doc/index.html) instead.

## `moodycamel::ConcurrentQueue<T>` by Cameron Desrochers

Available from <https://github.com/cameron314/concurrentqueue>, commit [`c680721`](https://github.com/cameron314/concurrentqueue/commit/c68072129c8a5b4025122ca5a0c82ab14b30cb03).

Affected files:

- `/deps/concurrentqueue/**`

Dual-licensed under Simplified BSD License and Boost Software License.

**NOTE** Added in [1.0.0](#v-1.0.0-section).

## `BS::thread_pool` by Barak Shoshany

Available from <https://github.com/bshoshany/thread-pool>, commit [`bd4533f`](https://github.com/bshoshany/thread-pool/commit/bd4533f1f70c2b975cbd5769a60d8eaaea1d2233) at tag [`v5.1.0`](https://github.com/bshoshany/thread-pool/releases/tag/v5.1.0).
  
Affected files:

- `/deps/thread-pool/**`
  
Under the MIT License.

**NOTE** Added in [1.1.0](#v-1.1.0-section).

## `{fmt}` by Victor Zverovich

Available from <https://github.com/fmtlib/fmt>, release [`12.1.0`](https://github.com/fmtlib/fmt/releases/tag/12.1.0).

Affected files:

- `/deps/slim_fmt/**`

Under an MIT-like License.

**NOTE** Added in [1.1.1](#v-1.1.1-section).

## PCG, Minimal C Implementation by Melissa O'Neill

Available from <https://www.pcg-random.org/download.html>. Currently developed at <https://github.com/imneme/pcg-c-basic>.

Affected files:

- `src/libam_support/ds/pcg_32_c.hh`

Licensed under the Apache 2.0 License.

**NOTE** Added in [1.3.0](#v-1.3.0-section).

## GNU Scientific Library by M. Galassi _et al._

Available from <https://www.gnu.org/software/gsl/>. The code used were adapted from the C version of [GSL 2.8 source tarball](https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz).

Affected files:

- `src/libam_support/ds/GslDiscreteDistribution.hh`, which is a re-implemented in C++ from [`randist/discrete.c`](https://cgit.git.savannah.gnu.org/cgit/gsl.git/tree/randist/discrete.c) in GSL 2.8.

Licensed under the GPL 3.0 license.

**NOTE** Added in [1.1.2](#v-1.1.2-section).

## qreverse by Wunkolo

Available from <https://github.com/Wunkolo/qreverse/>, commit [`fd54e4b`](https://github.com/Wunkolo/qreverse/commit/fd54e4b59e2e0b5051ffcfdab9878f10815bb8cd).

Affected files:

- `src/libam_support/seq/qreverse.c`
- `src/libam_support/seq/qreverse.h`

Licensed under the MIT license.

**NOTE** Added in [1.3.0](#v-1.3.0-section).
