(art_profile_builder-usage-section)=
# Usage of `art_profile_builder`

This new executable is designed to suppress the old `art_profiler_illumina` Shell/Perl scripts for building ART-compatible quality profiles. It supports building profiles from FASTQ and SAM/BAM files. Added in [1.1.9](#v-1.1.9-section).

**NOTE** All reads with quality will be used to build the profile. So if you use SAM/BAM files as input, please make sure that primiary- and un-aligned reads have quality information, and secondary- or supplementary-aligned reads have no quality information.

**NOTE** For MPI users: Currently, this executable runs under one slot (`mpiexec -n 1`) only. If more than one slot is provided, processes with MPI rank larger than 0 will terminate themselves.

## Options

- `--i-file`: Input file path. FASTQ/SAM/BAM files are supported. We currently rely on HTSLib to determine file format automatically. See also: [`htsfile(1)`](https://www.htslib.org/doc/htsfile.html).
- `--i-num_threads`: Number of threads used to decompress BAM input.
- `--read_len`: Read length. For each read, only the first `read_len` bases are used. If the read is shorter than `read_len`, only the available bases are used.
- `--is_pe`: For SAM/BAM input, whether the input is paired-end.
- `--o-file1`: Output quality profile file path for read 1.
- `--o-file2`: Output quality profile file path for read 2. Required if `--is_pe` is set.
- `--parallel`: Same as [`art_modern`'s `--parallel`](#parallelism-section).
- `--old_behavior`: Simulate the behavior of original ART profile builder. If set, all qualities will be offset by 1. See also: [The Bug in `Srcurity.md`](#original-art-profile-builder-bug).

## Environment Variables

Same as [`art_modern`'s environment variables](#am-environment-variables-section).

## Example

Building profiles from single-end FASTQ files:

```shell
art_profile_builder \
    --i-file out_se.fq \
    --read_len 36 \
    --o-file1 out_se_art_cxx_fq.txt \
    --parallel 2 \
    --i-num_threads 4
```

Building profiles from paired-end FASTQ files:

```shell
for i in 1 2; do
    art_profile_builder \
        --i-file out_pe."${i}".fq \
        --read_len 36 \
        --o-file1 out_pe_art_cxx_fq_R"${i}".txt \
        --parallel 2 \
        --i-num_threads 4
done
```

Building profiles from single-end SAM/BAM files:

```shell
art_profile_builder \
    --i-file out_se.sam \
    --read_len 36 \
    --o-file1 out_se_art_cxx_sam.txt \
    --parallel 2 \
    --i-num_threads 4
```

Building profiles from paired-end SAM/BAM files:

```shell
art_profile_builder \
    --i-file out_pe.sam \
    --read_len 36 \
    --is_pe \
    --o-file1 out_pe_art_cxx_sam_R1.txt \
    --o-file2 out_pe_art_cxx_sam_R2.txt \
    --parallel 2 \
    --i-num_threads 4
```
