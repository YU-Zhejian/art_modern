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

Available from <https://github.com/abseil/abseil-cpp>, release [`20250814.1`](https://github.com/abseil/abseil-cpp/releases/tag/20250814.1). The bundled source code is tailored from the [release tarball](https://github.com/abseil/abseil-cpp/releases/download/20250814.0/abseil-cpp-20250814.0.tar.gz).

Affected files:

- `/deps/slim_abseil/**`

Under the Apache-2.0 license.

**NOTE** Added in 1.1.0.

## `{fmt}` by Victor Zverovich

Available from <https://github.com/fmtlib/fmt>, release [`12.0.0`](https://github.com/fmtlib/fmt/releases/tag/12.0.0).

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

## VMT19937 by Fabio Cannizzo

Available from <https://github.com/fabiocannizzo/VMT19937>, commit [fa93111](https://github.com/fabiocannizzo/VMT19937/commit/fa93111bfc6f56f25c315990430a3487cdff9935).

Affected files:

- `/deps/VMT19937-master/**`

Licensed under the MIT license.

**NOTE** Added in 1.1.8.

## `xoshiro` by Nessan Fitzmaurice, David Blackman and Sebastiano Vigna

Available from <https://github.com/nessan/xoshiro> under tag [v1.0.0](https://github.com/nessan/xoshiro/releases/tag/v1.1.0) and commit [176fa19](https://github.com/nessan/xoshiro/commit/176fa191c8493e4c5cb06a44bc083010664fe39b).

Affected files:

- `/deps/xoshiro/**`

The version we used (`vigna.h`) is a C++ implementation migrated from the C code of original authors (David Blackman and Sebastiano Vigna) as we're currently not supporting C++20. The original C code is distributed under the [CC-0 license](http://creativecommons.org/publicdomain/zero/1.0/).

Other files are licensed under the MIT license.

**NOTE** Added in 1.1.8.

## Other Random Generators Recommended and Adapted by M.E. O'Neill

M.E. O'Neill, the author of the PCG random generator, had recommended several other random generators in his [blog post](https://www.pcg-random.org/posts/some-prng-implementations.html). He generously adapted those random generators from original authors and made them available under the MIT license.

Affected files in `deps/other_rngs` with their original GitHub Gist URLs:

- [`jsf.hpp`](https://gist.github.com/imneme/85cff47d4bad8de6bdeb671f9c76c814)
- [`gjrand.hpp`](https://gist.github.com/imneme/7a783e20f71259cc13e219829bcea4ac)
- [`sfc.hpp`](https://gist.github.com/imneme/f1f7821f07cf76504a97f6537c818083)
- [`lehmer.hpp`](https://gist.github.com/imneme/aeae7628565f15fb3fef54be8533e39c)
- [`splitmix.hpp`](https://gist.github.com/imneme/6179748664e88ef3c34860f44309fc71)
- [`arc4.hpp`](https://gist.github.com/imneme/4f2bf4b4f3a221ef051cf108d6b64d5a)

Licensed under the MIT license.

**NOTE** Added in 1.1.8.
