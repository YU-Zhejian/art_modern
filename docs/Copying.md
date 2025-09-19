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

**NOTE** This project uses the source code retrieved from [GitHub 1.22.1 Release Tarball](https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2).

## `gitignore` by GitHub

Available from <https://github.com/github/gitignore>.

Affected files:

- `/.gitignore`

With CC0 1.0 Universal License.

**NOTE** Added in 1.0.0.

## `libceu` by YU Zhejian

Affected files:

- `/deps/cmake_collections/**`
- `/deps/slim_libceu/**`
  
With MIT License.

**NOTE** Added in 1.0.0.

**NOTE** This project was dead. New projects should consider [Boost.Predef](https://www.boost.org/doc/libs/1_87_0/libs/predef/doc/index.html) instead.

## `moodycamel::ConcurrentQueue<T>` by Cameron Desrochers

Available from <https://github.com/cameron314/concurrentqueue>, commit [`c680721`](https://github.com/cameron314/concurrentqueue/commit/c68072129c8a5b4025122ca5a0c82ab14b30cb03).

Affected files:

- `/deps/concurrentqueue/**`

Dual-licensed under Simplified BSD License and Boost Software License.

**NOTE** Added in 1.0.0.

## `BS::thread_pool` by Barak Shoshany

Available from <https://github.com/bshoshany/thread-pool>, commit [`aa3fbfb`](https://github.com/bshoshany/thread-pool/commit/aa3fbfbe80762fe3ac90e2bf05e153b92536277a) at tag [`v5.0.0`](https://github.com/bshoshany/thread-pool/releases/tag/v5.0.0).
  
Affected files:

- `/deps/thread-pool/**`
  
Under the MIT License.

**NOTE** Added in 1.1.0.

## Abseil by Abseil Authors

Available from <https://github.com/abseil/abseil-cpp>, release [`20250814.0`](https://github.com/abseil/abseil-cpp/releases/tag/20250814.0). The bundled source code is tailored from the [release tarball](https://github.com/abseil/abseil-cpp/releases/download/20250814.0/abseil-cpp-20250814.0.tar.gz).

Affected files:

- `/deps/slim_abseil/**`

Under the Apache-2.0 license.

**NOTE** Added in 1.1.0.

## `{fmt}` by Victor Zverovich

Available from <https://github.com/fmtlib/fmt>, commit [`40626af`](https://github.com/fmtlib/fmt/commit/40626af88bd7df9a5fb80be7b25ac85b122d6c21) (tag [`11.2.0`](https://github.com/fmtlib/fmt/releases/tag/11.2.0)).

Affected files:

- `/deps/slim_fmt/**`

Under an MIT-like License.

**NOTE** Added in 1.1.1.

## PCG by Melissa O'Neill

Available from <https://www.pcg-random.org/downloads/pcg-cpp-0.98.zip>. Currently developed at <https://github.com/imneme/pcg-cpp>.

Affected files:

- `/deps/pcg-cpp-0.98/**`

Licensed under the Apache 2.0 License.

**NOTE** Added in 1.1.1.

## GNU Science Library by M. Galassi _et al._

Available from <https://www.gnu.org/software/gsl/>. The code used were adapted from the C version of [GSL 2.8 source tarball](https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz).

Affected files:

- `src/libam_support/ds/GslDiscreteDistribution.hh`, which is a re-implemented in C++ from <randist/discrete.c> in GSL 2.8.

Licensed under the GPL 3.0 license.

**NOTE** Added in 1.1.2.
