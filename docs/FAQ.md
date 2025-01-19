# Frequently Asked Questions (FAQs)

## Simulation

### What kinds of CPUs are supported by the simulator?

Although this application should theoretically support both endianness, only little endian is tested. That is, if you're working on an Intel or AMD CPU, this application should work fine.

### How to split produced pair-end/mate-pair sequencing results to 2 FASTQ files?

This can be done through [`seqtk`](https://github.com/lh3/seqtk). For example, to split `1.fq`:

```shell
# Read 1
seqtk seq 1.fq -1 > 1_1.fq
# Read 2
seqtk seq 1.fq -2 > 1_2.fq
```

You may also generate SAM/BAM files and extract PE FASTQ from them using `samtools`.

### Is there an interface for Python, Java, and more?

We may develop a C interface in the future after the APIs and design of the core library are settled.

### How to add adaptors \& primers to the reads?

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adapters to them, and use `templ` mode to introduce sequencing errors.

## Software Engineering

### Why you use `boost::filesystem` instead of `std::filesystem`, which was introduced in C++17?

The `std::filesystem` implementation in different compilers is not consistent. Some may require additional linker flags (`-lstdc++fs` for GCC <= 9.1; `-lc++fs` for Clang <= 9.0; `-lc++experimental` for Clang <= 7.0). There are also various reports in how those implementations deal with the terminating `/` when invoking `std::filesystem::creare_directories()`. So for the sake of simplicity, we use `boost::filesystem` instead.

See [this note in `cppreference`](https://en.cppreference.com/w/cpp/filesystem), [this StackOverflow question](https://stackoverflow.com/questions/53365538/how-to-determine-whether-to-use-filesystem-or-experimental-filesystem) and [this AskUbuntu question](https://askubuntu.com/questions/1256440/how-to-get-libstdc-with-c17-filesystem-headers-on-ubuntu-18-bionic).

### Does this program support cross-compilation?

Currently, it does not due to the extensive use of `test_run` in CMake scripts.

### I am using platforms other than x86\_64 (e.g., ARM, RISC-V, LoongArch, etc.), and spotted a bug!

Please submit a bug-report with instructions on how I may emulate your platform using free emulators like [QEMU](https://www.qemu.org/).

## My question is still unanswered

Please feel free to open an issue on GitHub or E-mail the maintainer.
