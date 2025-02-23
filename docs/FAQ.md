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

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adaptors to them, and use `templ` mode to introduce sequencing errors.

## Software Engineering

### Why you use `boost::filesystem` instead of `std::filesystem`, which was introduced in C++17?

The `std::filesystem` implementation in different compilers is not consistent. Some may require additional linker flags (`-lstdc++fs` for GCC <= 9.1; `-lc++fs` for Clang <= 9.0; `-lc++experimental` for Clang <= 7.0). There are also various reports in how those implementations deal with the terminating `/` when invoking `std::filesystem::creare_directories()`. So for the sake of simplicity, we use `boost::filesystem` instead.

See [this note in `cppreference`](https://en.cppreference.com/w/cpp/filesystem), [this StackOverflow question](https://stackoverflow.com/questions/53365538/how-to-determine-whether-to-use-filesystem-or-experimental-filesystem) and [this AskUbuntu question](https://askubuntu.com/questions/1256440/how-to-get-libstdc-with-c17-filesystem-headers-on-ubuntu-18-bionic).

### Why don't you use Intel Thread-Building Blocks (TBB)?

Intel TBB is an excellent library that supports diverse parallel programming models and data structures. However, this library loads into the memory in run-time, which would consume a lot of time in short-running applications. This also makes the project incapable of being distributed in a fully static form (See [Design.md](Design.md)).

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
  - If you use Intel (R) oneAPI DPC++/C++ Compiler, it will be `icpx --version`. You may need to set environment variables through e.g., `source /opt/intel/oneapi/setvars.sh`, before invoking this command.
  - Consult your compiler's documentation for other compilers.

## My question is still unanswered

Please feel free to open an issue on GitHub or E-mail the maintainer.
