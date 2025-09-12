# Copyright

This project uses GPL v3 License, a copy of which is available at [License.md](../License.md).

---

This project is based on the source code of [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) by [Weichun Huang](mailto:whduke@gmail.com) _et al._, under [GNU GPL v3](https://www.gnu.org/licenses/) license

---

This project uses code from the following projects (in no particular order):

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

Available from <https://github.com/YU-Zhejian/libceu>, commit `d1c043d8`.

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

- `/deps/thread_pool/**`
  
Under the MIT License.

**NOTE** Added in 1.1.0.

## Abseil by Google Inc.

Available from <https://github.com/abseil/abseil-cpp>, commit `fd8b35b9`.

Affected files:

- `/deps/slim_abseil/**`

**NOTE** Added in 1.1.0.

## Bravo C++ B-tree by Google Inc. \& German Mendez

Available from <https://github.com/Kronuz/cpp-btree>, commit `3e9f417`.

Affected files:

- `/deps/cpp-btree/**`

Under the Apache License, Version 2.0.

**NOTE** Added in 1.1.1.

**NOTE** No longer used in main code since 1.1.2.

## `{fmt}` by Victor Zverovich

Available from <https://github.com/fmtlib/fmt>, commit `a3d05d70` (tag `7.1.3`).

Affected files:

- `/deps/slim_libfmt/**`

Under an MIT-like License.

**NOTE** Added in 1.1.1.

**NOTE** This is the last version that does not come with C++20 modules.

## PCG by Melissa O'Neill \& PCG Project Contributors

Available from <https://github.com/imneme/pcg-cpp>.

Affected files:

- `/deps/pcg-cpp-0.98/**`

Licensed under the MIT License and Apache 2.0 License.

**NOTE** Added in 1.1.1.

**NOTE** This project used the version retried [here](https://www.pcg-random.org/downloads/pcg-cpp-0.98.zip), which was at that time licensed under Apache 2.0 License.

## `cgranges` by Dana-Farber Cancer Institute

Available from <https://github.com/lh3/cgranges>, commit `b3d5e2c5`.

Affected files:

- `deps/cgranges/*`.

With MIT License.

**NOTE** Added in 1.1.5.

## `EWAHBoolArray` by Daniel Lemire

Available from <https://github.com/lemire/EWAHBoolArray>, commit `a6a28d2`.

Affected files:

- `deps/EWAHBoolArray/*`.

With Apache/MIT License.

**NOTE** Added in 1.1.5.
