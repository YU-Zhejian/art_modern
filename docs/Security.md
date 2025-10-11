# Security Policy

## Reporting a Vulnerability

Send an E-mail to maintainer in private if vulnerabilities were found.

## Known Bugs in `art_modern`

% TODO

## Known Bugs in Original ART

### Original ART may produce illegal SAM files.

### Original ART Profile Builder Problem

Steps to reproduce:

Given file `noqual_test/1.1.fq`:

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

Create a profile with command:

```shell
deps/ART_profiler_illumina/art_profiler_illumina \
    opt/noqual_test/noqual_test_ \
    opt/noqual_test \
    fq
```

```shell
printf '>a\nAGCTAGCTACAGCTAGCTACAGCTAGCTACAGCTAGCTACAGCTAGCTACAGCTAGCTACTGATCG\n' > opt/noqual_test/ref.fa
art_illumina \
    --in opt/noqual_test/ref.fa \
    --out opt/noqual_test/art_out \
    --qprof1 opt/noqual_test/noqual_test_R1.txt \
    --rcount 100 --len 20 -amp --samout --maskN 0
```

The qualities of generated file will be `"`, instead of `!` in the input.
