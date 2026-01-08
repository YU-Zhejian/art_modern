# Frequently Asked Questions (FAQs)

## Simulation

(fastq-split-section)=
### How to split produced pair-end/mate-pair sequencing results to 2 FASTQ files?

This can be done through [`seqtk`](https://github.com/lh3/seqtk). For example, to split `1.fq`:

```shell
# Read 1
seqtk seq a.fq -1 > a_1.fq
# Read 2
seqtk seq a.fq -2 > a_2.fq
```

You may also generate SAM/BAM files and extract PE FASTQ from them using `samtools`. For example:

```shell
samtools fastq -1 a_1.fq -2 a_2.fq -N a.sam
```

**NOTE** Please sort the generated FASTQ files using `seqkit sort` to ensure reads with the same name are in the same order.

### How to add adaptors \& primers to the reads?

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adaptors to them, and use template mode to introduce sequencing errors.

### How to support new Illumina models?

See the documentation of [`art_profile_builder`](#art_profile_builder-usage-section) to build new profiles from real sequencing data.

## Software Engineering

### What kinds of CPUs are supported by the simulator?

Although this application should theoretically support both endianness, only little endian is tested. That is, if you're working on an Intel or AMD CPU, this application should work fine.

### Does this program support cross-compilation?

Currently, it does not due to the extensive use of `test_run` in CMake scripts.

### I am using platforms other than x86\_64 (e.g., ARM, RISC-V, LoongArch, etc.), and spotted a bug

Please submit a bug-report with instructions on how I may emulate your platform using free emulators like [QEMU](https://www.qemu.org/).

## Community

### I am new to GNU/Linux. How can I gather required information for a bug report?

Where you may find the required information:

- Name and version of your current operating system in a human-readable way.

  You may find that out by reading `/etc/lsb-release` or use tools like [neofetch](https://github.com/dylanaraps/neofetch), [screenFetch](https://github.com/KittyKatt/screenFetch), or [fastfetch](https://github.com/fastfetch-cli/fastfetch).

- Version of the kernel.

  Easy, just run `uname -a`, or `cat /proc/version`.

- Version of the compiler.

  - If your CMake correctly identifies your compiler, you may put the first several lines of CMake output to the bug report.
  - If you use GCC, it will be `g++ --version`.
  - If you use LLVM Clang, it will be `clang++ --version`.
  - If you use Intel (R) OneAPI DPC++/C++ Compiler, it will be `icpx --version`. You may need to set environment variables through e.g., `source /opt/intel/oneapi/setvars.sh` (if you install it in default location using APT), before invoking this command.
  - Consult your compiler's documentation for other compilers.

## My question is still unanswered

Please feel free to open an issue on GitHub or E-mail the maintainer.
