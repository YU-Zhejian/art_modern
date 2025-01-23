---
name: Bug report
about: Create a report to help improve this package
title: ''
labels: ''
assignees: ''

---

## Session Information

Please provide the following information about your system by replacing the placeholders. Avoid any blanks that may expose your personal information by deleting the placeholder while unticking the box. You may use `/workbench` as a placeholder to your working directory.

- [ ] Name and version of your current operating system in an human-readable way.

  `Linux Mint 22 Wilma`

- [ ] Version of the kernel.

  ```text
  Linux SOME-WSL 5.15.167.4-microsoft-standard-WSL2 #1 SMP Tue Nov 5 00:21:55 UTC 2024 x86_64 x86_64 x86_64 GNU/Linux
  ```

- [ ] Version of the compiler.

  ```text
  gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0
  ```

- [ ] Version of dependencies and how they're installed. List those in the following format while removing those unused ones.

  - GNU Binutils; `2.42-4ubuntu2.3`; Official APT Source.
  - GNU Science Library; `2.7.1==he838d99_0`; Conda under `conda-forge` channel.
  - Boost; `1.87.0`; Compiled from official source tarballs.
  - HTSLib; `1.21`; Bundled.
  - Intel Math Kernel Library; `2023.0.0`; From official binary distribution.
  - LLVM; `14.2.0`; Compiled from source against musl libc with instructions from [here](https://wiki.musl-libc.org/building-llvm.html).
  - GNU Compiler Collection; `14.2.0`; Unknown origin.

- [ ] `art_modern --version` output.

  ```text
  [2024-12-02 00:16:21.329146] info: YuZJ Modified ART_Illumina (art_modern v. 1.0.0)
  [...]
  ART: 2.5.8, ART_MODERN: 1.0.0
  ART_MODERN_LINK_LIBS: Boost::filesystem;[...];libceu
  USING HTSLib: labw_slim_htslib, ver. 1.15
  [...]
            bool        (1 size):      0 -> 1
  [2024-12-02 00:16:21.360571] info: EXIT
  ```

- [ ] The CMake commandline you used to build the code.

  ```shell
  env -C opt/build_release cmake -DCMAKE_BUILD_TYPE=Release -G Ninja "$(pwd)"
  env -C opt/build_release ninja -j40
  ```

- [ ] The CMake configuration log, if possible.

  ```text
  -- The C compiler identification is GNU 14.2.1
  -- The CXX compiler identification is GNU 14.2.1
  -- Detecting C compiler ABI info
  -- Detecting C compiler ABI info - done
  [...]
  -- Looking for sys/utsname.h - found
  -- Looking for _mingw.h
  -- Looking for _mingw.h - not found
  ART_MODERN_LINK_LIBS=Boost::filesystem;[...];libceu
  Testing disabled since CEU_CM_SHOULD_ENABLE_TEST was set to FASE
  ```

- [ ] The compilation log, if possible.

  ```text
  [86/86] Linking CXX executable art_modern
  ```

- [ ] The actual libraries linked to the executable. You may find that out using `ldd` or `readelf -d`.

  ```text
  linux-vdso.so.1 (0x00007ffc101e0000)
  libart_modern_lib.so => /workbench/opt/build_release/libart_modern_lib.so (0x00007f1cd746d000)
  libboost_regex.so.1.85.0 => /lib64/libboost_regex.so.1.85.0 (0x00007f1cd741f000)
  [...]
  libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00007f1cd6bab000)
  libc.so.6 => /lib64/libc.so.6 (0x00007f1cd6999000)
  /lib64/ld-linux-x86-64.so.2 (0x00007f1cd750d000)
  ```

  ```text
  Dynamic section at offset 0x69c88 contains 47 entries:
    Tag        Type                         Name/Value
  0x0000000000000001 (NEEDED)             Shared library: [libart_modern_lib.so]
  0x0000000000000001 (NEEDED)             Shared library: [libboost_regex.so.1.85.0]
  0x0000000000000001 (NEEDED)             Shared library: [libboost_program_options.so.1.85.0]
  [...]
  0x0000000000000001 (NEEDED)             Shared library: [libgcc_s.so.1]
  0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
  0x000000000000001d (RUNPATH)            Library runpath: [/workbench/opt/build_release:/workbench/opt/build_release/deps/labw_slim_htslib:/workbench/opt/build_release/deps/slim_libceu]
  ```

- [ ] The git commit hash of the code you used.

  ```text
  d24316255d794afb6c0334686b261489f4a7d9c7
  ```

**NOTE: IF YOU'RE IN A CONTROLLED ACCESS ENVIRONMENT PLEASE MAKE SURE THE INFORMATION YOU SEND DOES NOT CONTAIN ANYTHING THAT MAY BE CONSIDERED CONFIDENTIAL!**

**NOTE: DO NOT SEND THE CMAKE BUILD FORDER AS IT MAY CONTAIN SENSITIVE INFORMATION!**

## Optional Questions

- All test cases can be passed (Set `-DCEU_CM_SHOULD_ENABLE_TEST=ON`, and run tests using `ctest`).
  - [ ] Yes
  - [ ] No
  - [ ] I don't know
- `make testsmall` can be passed.
  - [ ] Yes
  - [ ] No
  - [ ] I don't know
- The bug can be reproduced using the most recent GCC available in your system under `Debug` preset of CMake (`-DCMAKE_BUILD_TYPE=Debug`), non-native code (`-DCEU_CM_SHOULD_USE_NATIVE=OFF`) and `STL` random generator (`-DUSE_RANDOM_GENERATOR=STL`).
  - [ ] Yes
  - [ ] No
  - [ ] N/A (E.g., problem arose in tests)
  - [ ] I don't know
- The bug can be reproduced without parallelization (`-DUSE_THREAD_PARALLEL=NOP` CMake option).
  - [ ] Yes
  - [ ] No
  - [ ] N/A (E.g., problem arose in tests)
  - [ ] I don't know

## Describe the Bug Below

### Expected Behavior

### Actual Behavior

### Minimal Input and Steps to Reproduce

(Attaching links to external files or script that transforms a file from an external source is highly encouraged)

### Hypothesized Cause
