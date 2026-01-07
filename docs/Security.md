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
- External HTSLib and NCBI NGS library incompatibilities: Both library defines `string_dup`.

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
