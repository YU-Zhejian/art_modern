# Copyright

This project uses GPL v3 License, a copy of which is available at [License.md](../License.md).

---

This project is based on the source code of [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) by [Weichun Huang](mailto:whduke@gmail.com) _et al._, under [GNU GPL v3](https://www.gnu.org/licenses/) license

---

This project uses code from the following projects (in no particular order). All code used is either in the public domain or under a license [compatible with GPL v3](https://www.gnu.org/licenses/license-list.en.html#GPLCompatibleLicenses).

## HTSlib by Genome Research Ltd.

Available from <https://github.com/samtools/htslib>.

Affected files:

- `/deps/labw_slim_htslib/cram/**`, under the MIT/Expat License.
- `/deps/labw_slim_htslib/**`, under The Modified-BSD License.

**NOTE** This project uses the source code retrieved from [SourceForge 1.21 Release Tarball](https://sourceforge.net/projects/samtools/files/samtools/1.21/htslib-1.21.tar.bz2/download).

## `gitignore` by GitHub

Available from <https://github.com/github/gitignore>

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

Available from <https://github.com/cameron314/concurrentqueue>, commit `6dd38b8`.

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

## Abseil by Google Inc.

Available from <https://github.com/abseil/abseil-cpp>, commit `fd8b35b9`.

Affected files:

- `/deps/slim_abseil/**`

Under the Apache-2.0 license.

**NOTE** Added in 1.1.0.

## `{fmt}` by Victor Zverovich

Available from <https://github.com/fmtlib/fmt>, commit [`7bdf062`](https://github.com/fmtlib/fmt/commit/7bdf0628b1276379886c7f6dda2cef2b3b374f0b) (tag [`7.1.3`](https://github.com/fmtlib/fmt/releases/tag/7.1.3)).

Affected files:

- `/deps/slim_fmt/**`

Under an MIT-like License.

**NOTE** Added in 1.1.1.

**NOTE** This is the last version that does not come with C++20 modules.

## PCG by Melissa O'Neill

Available from <https://www.pcg-random.org/downloads/pcg-cpp-0.98.zip>. Currently developed at <https://github.com/imneme/pcg-cpp>.

Affected files:

- `/deps/pcg-cpp-0.98/**`

Licensed under the Apache 2.0 License.

**NOTE** Added in 1.1.1.

## GNU Science Library by M. Galassi et al.

Available from <https://www.gnu.org/software/gsl/>. The code used were adapted from the C version of GSL 2.8 source tarball.

Affected files:

- `src/libam_support/ds/GslDiscreteDistribution.hh`, which is a re-implemented in C++ from <randist/discrete.c> in GSL 2.8.

Licensed under the GPL 3.0 license.

**NOTE** Added in 1.1.2.
