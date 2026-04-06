# Frequently Asked Questions (FAQs)

## Simulation

(fastq-split-section)=
### How to split the produced pair-end/mate-pair sequencing results into 2 FASTQ files?

This can be done through [`seqtk`](https://github.com/lh3/seqtk). For example, to split `a.fq`:

```shell
# Read 1
seqtk seq a.fq -1 > a_1.fq
# Read 2
seqtk seq a.fq -2 > a_2.fq
```

You may also generate SAM/BAM files and extract PE FASTQ from them using `samtools`. For example:

```shell
samtools fastq -1 a_1.fq -2 a_2.fq -N a.sam
seqkit sort a_1.fq > a_sorted_1.fq
seqkit sort a_2.fq > a_sorted_2.fq
```

**NOTE:** If you generate FASTQs through SAM, please sort the generated FASTQ files using `seqkit sort` to ensure reads with the same name are in the same order.

### How to add adaptors \& primers to the reads?

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adaptors, and use the `template` mode to introduce sequencing errors.

### How to support new Illumina models?

See the documentation of [`art_profile_builder`](#art_profile_builder-usage-section) to build new profiles from real sequencing data.

## Software Engineering

### What kinds of CPUs are supported by the simulator?

Although this application should theoretically support both endianness, only little endian is tested. That is, if you're working on an Intel or AMD CPU, this application should work fine.

### Does this program support cross-compilation?

Currently, it does not work due to the extensive use of `test_run` in CMake scripts.

### I am using platforms other than x86\_64 (e.g., ARM, RISC-V, LoongArch, etc.), and spotted a bug

Please submit a bug report with instructions on how I can emulate your platform using free emulators like [QEMU](https://www.qemu.org/).

## My question is still unanswered

Please feel free to open an issue on GitHub or E-mail the maintainer.
