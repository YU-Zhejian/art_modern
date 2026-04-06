# Security Policy

## Reporting a Vulnerability

Send a private email to the maintainer if vulnerabilities were found.

## Known Bugs in `art_modern`

See [News.md](News.md) for details of the bugs fixed in each release.

### Known Build-Time Incompatibilities

- Clang earlier than 9 with Boost later than 1.78 may raise bugs when using headers from `boost/math`, since Boost may misidentify Clang as an older GCC that does not support C++11. See [here](https://www.boost.org/doc/libs/1_87_0/libs/math/doc/html/math_toolkit/history2.html#math_toolkit.history2.math_3_0_0_boost_1_76) for more details.
- Boost versions earlier than 1.65 do not contain `boost/asio/thread_pool.hpp` or `boost/asio/post.hpp`. See [here](https://www.boost.org/doc/libs/1_66_0/doc/html/boost_asio/reference/thread_pool.html) for its first introduction. Under these circumstances, try `-DUSE_THREAD_PARALLEL=BS` in CMake options (introduced below) as an alternative.
- GCC 13.3.0 on Haiku OS hrev58590 may cause a kernel panic that locks up the entire system during builds.
- GCC would fail on Debian GNU/Hurd.
- External HTSLib and NCBI NGS library incompatibilities: Both library defines `string_dup`.
- Builds that use Intel MKL or Intel C++/DPC++ compilers may be incompatible with [Citrix Workspace app](https://www.citrix.com/downloads/workspace-app/) 25.08.10.111 under Linux.

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

Create a profile with the command:

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

The qualities of the generated file will be `"` (Phred-score 1), instead of `!` (Phred-score 0) in the files that generate the profile.
