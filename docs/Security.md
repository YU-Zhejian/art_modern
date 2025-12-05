# Security Policy

## Reporting a Vulnerability

Send an E-mail to maintainer in private if vulnerabilities were found.

## Known Bugs in `art_modern`

See [News.md](News.md) for details of the bugs fixed in each release.

### Known Build-Time Incompatibilities

- Clang earlier than 9 with Boost later than 1.78 may raise bug when using headers from `boost/math` as boost may misidentify Clang as GCC that does not support C++11. See [here](https://www.boost.org/doc/libs/1_87_0/libs/math/doc/html/math_toolkit/history2.html#math_toolkit.history2.math_3_0_0_boost_1_76) for more details.
- Boost earlier than 1.65 does not contain `boost/asio/thread_pool.hpp` or `boost/asio/post.hpp`. See [here](https://www.boost.org/doc/libs/1_66_0/doc/html/boost_asio/reference/thread_pool.html) for its first introduction. Under this circumstance, try `-DUSE_THREAD_PARALLEL=BS` in CMake options (introduced below) as an alternative.
- GCC 13.3.0 on Haiku OS hrev58590 may generate a kernel panic that jam the entire system while building.
- GCC would fail on Debian GNU/Hurd.
- CYGWIN users: Clang 20.1.8 and Boost 1.66.0 build failure on commit `c996d6071c076b50b3b73c364989b12e80685358`.

    ```text
    Compile-time C std.: ver. unknown ISO C (__STDC_VERSION__ undefined)
    Compile-time C++ std.: ver. 17 (201703)
    Compiled at: N/A due to reproducible build
    Compiler Identification:
            Clang compatible version number: 20.1.8 (20.1.8 (https://cygwin.com/git/cygwin-packages/clang e8422310037f1f8bc3ec6144c1d5ef19eadb28fc))
            GCC compatible version number: 4.2.1
            __VERSION__ string: Clang 20.1.8 (https://cygwin.com/git/cygwin-packages/clang e8422310037f1f8bc3ec6144c1d5ef19eadb28fc)
    Compile-time C Types max, min, etc. limits:
            char           (1 size):      -128 -> +127
            schar          (1 size):      -128 -> +127
            uchar          (1 size):      +0 -> 255
            size_t         (8 size):      +0 -> 18446744073709551615
            ptrdiff_t      (8 size):      -9223372036854775808 -> +9223372036854775807
            short          (2 size):      -32768 -> +32767
            ushort         (2 size):      +0 -> +65535
            int            (4 size):      -2147483648 -> +2147483647
            uint           (4 size):      +0 -> 4294967295
            long           (8 size):      -9223372036854775808 -> +9223372036854775807
            ulong          (8 size):      +0 -> 18446744073709551615
            llong          (8 size):      -9223372036854775808 -> +9223372036854775807
            ullong         (8 size):      +0 -> 18446744073709551615
            bool           (1 size):      +0 -> +1
    Compile-time OS info:
            PRIMARY OS='CYGWIN'
            CYGWIN API ver. 0.357, with dll (cygwin1) ver. 3006.5
            POSIX.1 Version: 202405
            POSIX.2 Version: 202405
            Single UNIX Specification (SUS) Version: 700
    Run-time OS info:
            POSIX UTSINFO:
                    sysname=CYGWIN_NT-10.0-26200
                    nodename=DESKTOP-CGV4O9H
                    release=3.6.5-1.x86_64
                    version=2025-10-09 17:21 UTC
                    machine=x86_64
            Windows OS Version Info:
                    Platform ID: Windows 11 25H2 (2025 Update)
                    Major Version: 10.0.26200
                    Service Pack: not installed
    ```

## Known Bugs in Original ART

### Original ART may Produce Illegal SAM Files

(original-art-profile-builder-bug)=
### Original ART Profile Builder Problem

The original ART profile builder would generate profiles with quality scores offset by 1.

Steps to reproduce:

Given file `noqual_test/1.1.fq` where qualities are all `!` (Phred-score 0):

```text
@1/1
AGCTAGCTAGCTAGCTAGCT
+
!!!!!!!!!!!!!!!!!!!!
```

and file `noqual_test/1.2.fq`:

```text
@1/2
AGCTAGCTAGCTAGCTAGCT
+
!!!!!!!!!!!!!!!!!!!!
```

and file `opt/noqual_test/ref.fa`:

```text
>a
AGCTAGCTACAGCTAGCTACAGCTAGCTACAGCTAGCTACAGCTAGCTACAGCTAGCTACTGATCG
```

Create a profile with command:

```shell
art_profiler_illumina opt/noqual_test/noqual_test_ opt/noqual_test fq
```

And then use the original ART to generate reads with the created profile:

```shell
art_illumina \
    --in opt/noqual_test/ref.fa \
    --out opt/noqual_test/art_out \
    --qprof1 opt/noqual_test/noqual_test_R1.txt \
    --rcount 2 --len 20 -amp --samout --maskN 0 --rndSeed 0
```

The generated file will look like:

```text
@a-2
TGCGTATTCGGCCTTTGTTT
+
""""""""""""""""""""
@a-1
CGTAATTATACGCGCATGAG
+
""""""""""""""""""""
```

The qualities of generated file will be `"` (Phred-score 1), instead of `!` (Phred-score 0) in the files where the profile is generated from.
