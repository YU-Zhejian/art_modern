# Copyright

This project uses GPL v3 License, a copy of which is available at [License.md](../License.md).

Except those acknowledged in Readme, this project uses code from the following projects:

## GitHub `gitignore` Project

Available from <https://github.com/github/gitignore>

Affected files:

- `/.gitignore`

With CC0 1.0 Universal License.

## YU Zhejian's `libceu` Project

Available from <https://github.com/YU-Zhejian/libceu>, commit `d1c043d8`.

Affected files:

- `/deps/cmake_collections/**`
- `/deps/slim_libceu/**`
  
With MIT License.

**NOTE** This program was dead. New projects should consider [Boost.Predef](https://www.boost.org/doc/libs/1_87_0/libs/predef/doc/index.html) instead.

## Cameron Desrochers's `moodycamel::ConcurrentQueue<T>` Project

Available from <https://github.com/cameron314/concurrentqueue>, commit `6dd38b8`.

Affected files:

- `/deps/concurrentqueue/**`

Dual-licensed under Simplified BSD License and Boost Software License.

## William Sherif's Nibble And A Half Project

Originally from <https://github.com/superwills/NibbleAndAHalf>.

Affected files:

- `/deps/base64/**`

Under an MIT-like License.

**NOTE** The version used in this project comes from [MMseqs2](https://github.com/soedinglab/MMseqs2), commit `b804fbe3`, which is copyrighted by The MMseqs2 Development Team under the MIT License. We use this version instead of the latest one since this version is clean and header-only.

**NOTE** Some modifications were done to supress clang-tidy warnings.

## Barak Shoshany's `BS::thread_pool` Project

Available from <https://github.com/bshoshany/thread-pool>, commit `aa3fbfb`.
  
Affected files:

- `/deps/thread_pool/**`
  
Under the MIT License.

## Google Inc. Abseil Project

Available from <https://github.com/abseil/abseil-cpp>, commit `fd8b35b9`.

Affected files:

- `/deps/slim_abseil/**`

## Google Inc. and German Mendez Bravo C++ B-tree Project

Available from <https://github.com/Kronuz/cpp-btree>, commit `3e9f417`.

Affected files:

- `/deps/cpp-btree/**`

Under the Apache License, Version 2.0.

## Victor Zverovich's `{fmt}` Project

Available from <https://github.com/fmtlib/fmt>, commit `a3d05d70` (tag `7.1.3`).

Affected files:

- `/deps/slim_libfmt/**`

Under an MIT-like License.

**NOTE** This is the last version that does not come with C++20 modules.
