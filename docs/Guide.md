# Guide on Miscellaneous Topics

This guide covers several miscellaneous topics that may be of interest to users and developers of this software.

## Guide on FASTA Format

FASTA format can be parsed by all parsers. However, please keep in mind that for `memory` parser, the following kinds of FASTA files are supported:

```text
>one_line_fasta
AAAAAAAAAAAAAAAAA
>multi_line_fasta
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
>fasta_with_empty_sequence
>fasta_with_empty_sequence_with_newlines



>fasta_with_spaces_in_name some description here
AAAAAAAAA
```

Note that `fasta_with_empty_sequence_with_newlines` is **NOT** supported by [PacBio Formats](https://pacbiofileformats.readthedocs.io/en/13.0/FASTA.html) or [NCBI GenBank FASTA Specification](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) and [NCBI GenBank Submission Guidelines](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#files).

The following kinds of FASTA files are **NOT** supported:

```text
> some_sequence_without_a_name
AAAAAAAA
```

Also, note that all characters other than `ACGTacgt` will be regarded as `N`. We do **NOT** support IUPAC codes.

**NOTE** For `htslib` parser, identical line lengths (except the last line) inside a contig is assumed. That means the following FASTA file is legal for `htslib` parser:

```text
>chr1
AAAAAAAAAAAA
AAAA
>chr2
AAAAAAAAAAAA
AAAAAAAAAAAA
AAAAAAAAAAAA
AA
>chr3
AA
```

But the following is not:

```text
>chr2
AAAA
AAAAAAAAAA
AAAAAAAA
AA
```

**NOTE** For read names, only characters before the first whitespace characters (space ` `, tabs `\t`, etc.) are read. That is, the FASTA file:

```text
>chr1 some attrs
AAAAAATTTTTT
>chr2 more attrs
AAAAAATTTTTT
```

Will be parsed into identical data structure with:

```text
>chr1
AAAAAATTTTTT
>chr2
AAAAAATTTTTT
```

Using empty file or `/dev/null` as input is allowed since 1.1.9. It will generate empty FASTA/FASTQ/PWA files as output. Remember to specify `memory` or `stream` as `--i-parser` `--i_type` and do **NOT** use SAM/BAM output writer in this case (as SAM/BAM output writer will think you're using streamed input and raise an exception).

## Performance Hint

When building `art_modern`, set [`USE_HTSLIB`](#use-htslib-section) to the latest HTSLib available on your system.  Please also make sure that your HTSLib has been compiled with `-O3 -mtune=native -march=native` and linked with [libdeflate](https://github.com/ebiggers/libdeflate). Set [`CMAKE_BUILD_TYPE`](#cmake-build-type-section) to `Release` or `RelWithDebInfo`, and [`USE_RANDOM_GENERATOR`](#use-random-generator-section) to `ONEMKL` on Intel/AMD machines or `PCG` on other machines.

When executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.

SAM/BAM output writers are memory- and time-consuming due to compression. If you don't need SAM/BAM output, please don't enable it.

## Platform-Specific Notes

### Debian, Ubuntu, Linux Mint, or Other Debian-Based Distributions

For minimal build, install the following packages using APT:

```shell
apt-get install -y \
    build-essential g++ binutils \
    libboost-all-dev zlib1g-dev \
    make python3 cmake sed grep coreutils
```

Then work using `cmake` as usual. CMake will use bundled source codes for missing dependencies.

For full-featured build (with external libraries and MPI support), install the following additional packages using APT:

```shell
apt-get install --no-install-recommends -y \
    liblzma-dev libbz2-dev libhts-dev \
    pkgconf \
    libfmt-dev \
    libconcurrentqueue-dev \
    libabsl-dev \
    openmpi-bin libopenmpi-dev
```

And use corresponding CMake options to disable using bundled libraries.

### Alpine Linux

For building static binaries, you may need to install the following packages using APK:

```shell
apk add g++ binutils \
    boost-dev boost1.84-static icu-static \
    zlib-dev zlib-static \
    make cmake python3 \
    coreutils sed grep \
    xz-static xz-dev \
    bzip2-static bzip2-dev \
    libdeflate-static libdeflate-dev
```

Here, `icu-static` is added to support Boost staic linking.

**NOTE** Please install the correct version of Boost static library.

**NOTE** `coreutils` is **MANDATORY** -- Those shipped with BusyBox will **NOT** work.

### Apple Mac OS X

Install Xcode command-line tools using:

```shell
xcode-select --install
```

Alternatively, you may get the latest version [here](https://developer.apple.com/download/all). An Apple account is required.

Download CMake from [here](https://cmake.org/download/). You should modify PATH in `~/.bashrc` or `~/.zshrc` to include the `bin` directory of the CMake installation. If you install the DMG image, it will commonly be located in `/Applications/CMake.app/Contents/bin/cmake`.

See the following section to build Boost from source.

Apple Mac OS X comes with zlib and its development headers pre-installed.

You may also set up dependencies using [Conda](https://docs.conda.io), [MacPorts](https://www.macports.org/) or [HomeBrew](https://brew.sh).

### Installing Boost from Source

If your system Boost library does not exist (e.g., on brand-new Apple MacOS X), is too old (e.g., older than 1.65.0) or ABI incompatible (e.g., compiled with GCC but you want to use Clang/LLVM), you may install Boost from source. Here is an example of installing Boost 1.89.0 using Clang/LLVM toolchain to `"${HOME}"/opt/boost-1.89.0-clang`:

```shell
# Assume we're using Boost 1.89.0
wget https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.gz
tar xvzf boost_1_89_0.tar.gz
cd boost_1_89_0
./bootstrap.sh # Build b2. There's no point to build b2 using Clang.
./b2 install \
    toolset=clang \
    cxxflags="-stdlib=libc++" \
    linkflags="-stdlib=libc++ -fuse-ld=lld" \
    --prefix="${HOME}"/opt/boost-1.89.0-clang
```

This should build all required Boost libraries and a majority of optional ones.

And then you may use CMake to build this project through:

```shell
mkdir -p build
cd build
# Set -DBoost_DIR accordingly.
# Older CMake may have different behaviour.
cmake .. -DBoost_DIR="${HOME}"/opt/boost-1.89.0-clang
```
