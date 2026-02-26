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

Using empty file or `/dev/null` as input is allowed since [1.1.9](#v-1.1.9-section). It will generate empty FASTA/FASTQ/PWA files as output. Remember to specify `memory` or `stream` as `--i-parser` `--i_type` and do **NOT** use SAM/BAM output writer in this case (as SAM/BAM output writer will think you're using streamed input and raise an exception).

## Performance Hint

When building `art_modern`, set [`USE_HTSLIB`](#use-htslib-section) to the latest HTSLib available on your system. Please also make sure your HTSLib has been compiled with `-O3 -mtune=native -march=native` and linked with [`libdeflate`](https://github.com/ebiggers/libdeflate). Set [`CMAKE_BUILD_TYPE`](#cmake-build-type-section) to `Release` or `RelWithDebInfo`, and [`USE_RANDOM_GENERATOR`](#use-random-generator-section) to `ONEMKL` on Intel/AMD machines or `PCG`/`SYSTEM_PCG` on other machines.

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

Here, `icu-static` is added to support Boost when performing static linking.

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

### MSYS2 MinGW64

See [MSYS2 Environments](https://www.msys2.org/docs/environments/) for the differences between MSYS2 environments.

Use `pacman` to install the following packages:

```text
autoconf-wrapper
automake-wrapper
base
bc
diffstat
dos2unix
man-db
mingw-w64-x86_64-autotools
mingw-w64-x86_64-binutils
mingw-w64-x86_64-boost
mingw-w64-x86_64-cmake
mingw-w64-x86_64-crt-git
mingw-w64-x86_64-fmt
mingw-w64-x86_64-gcc
mingw-w64-x86_64-gcc-libs
mingw-w64-x86_64-gdb
mingw-w64-x86_64-headers-git
mingw-w64-x86_64-htslib
mingw-w64-x86_64-libwinpthread
mingw-w64-x86_64-pkgconf
mingw-w64-x86_64-samtools
mingw-w64-x86_64-winpthreads
p7zip
pactoys
parallel
patch
patchutils
pcre
pcre2
pkgfile
procps-ng
psmisc
reflex
rlwrap
sqlite
wcd
xorg-util-macros
xorgproto
zip
zsh
```

% TODO

### MSYS2 Clang64

Use `pacman` to install the following packages in addition:

```text
mingw-w64-clang-x86_64-autotools
mingw-w64-clang-x86_64-boost
mingw-w64-clang-x86_64-boost-libs
mingw-w64-clang-x86_64-clang
mingw-w64-clang-x86_64-cmake
mingw-w64-clang-x86_64-fmt
mingw-w64-clang-x86_64-htslib
mingw-w64-clang-x86_64-python
mingw-w64-clang-x86_64-samtools
```

### Installing Boost from Source

If your system Boost library does not exist (e.g., on brand-new Apple Mac OS X), is too old (e.g., older than 1.65.0) or ABI incompatible (e.g., compiled with GCC, but you want to use Clang/LLVM), you may install Boost from source. Here is an example of installing Boost 1.89.0 using Clang/LLVM toolchain to `"${HOME}"/opt/boost-1.89.0-clang`:

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
    --prefix="${HOME}"/opt/boost-1.89.0-clang \
    --ignore-site-config \
    variant=release \
    threading=multi
```

This should build all required Boost libraries and most optional ones.

And then you may use CMake to build this project through:

```shell
mkdir -p build
cd build
# Set -DBoost_DIR accordingly.
# Older CMake may have different behaviour.
cmake .. -DBoost_DIR="${HOME}"/opt/boost-1.89.0-clang
```

### Installing `{fmt}` from Source

```shell
wget https://github.com/fmtlib/fmt/releases/download/12.0.0/fmt-12.0.0.zip
unzip fmt-12.0.0.zip
cd fmt-12.0.0
mkdir -p build
cd build
# Build shared library. Set to OFF to build static library.
cmake .. \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX="${HOME}/opt/fmt-12.0.0-clang"
cmake --build . -j "$(nproc)" --target fmt
cmake --install .
```

And then you may use CMake to build this project through:

```shell
mkdir -p build
cd build
PKG_CONFIG_PATH="${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/:${PKG_CONFIG_PATH:-}" \
    cmake ..
```

### Installing NCBI NGS SDK from Source

Install NCBI VDB first.

```shell
axel https://github.com/ncbi/ncbi-vdb/archive/refs/tags/3.3.0.zip
unzip ncbi-vdb-3.3.0.zip
cd ncbi-vdb-3.3.0
./configure --prefix="${HOME}"/opt/ncbi-vdb-3.3.0
make -j $(nproc) all install
make install
cd ..
```

Install NCBI NGS SDK as a part of NCBI SRA Toolkit.

```shell
axel https://github.com/ncbi/sra-tools/archive/refs/tags/3.3.0.zip
unzip sra-tools-3.3.0.zip
cd sra-tools-3.3.0
./configure \
    --prefix="${HOME}"/opt/sra-tools-3.3.0 \
    --with-ncbi-vdb-prefix="${HOME}"/opt/ncbi-vdb-3.3.0 
make clean install -j $(nproc) BUILD_TOOLS_LOADERS=ON
cd ..
cp -r \
    "${HOME}"/opt/ncbi-vdb-3.3.0/include/ \
    "${HOME}"/opt/sra-tools-3.3.0/bin/ncbi/schema
```
