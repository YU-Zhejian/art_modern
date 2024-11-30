# Readme for `art_modern`: Modernized Artificial Read Transcription (ART) for Simulation of Diverse Next-Generation Sequencing Reads

## Quick Start

Build the project using:

```shell
mkdir -p build_release
env -C build_release cmake -DCMAKE_BUILD_TYPE=Release -DCEU_CM_SHOULD_ENABLE_TEST=FALSE ..
env -C build_release make -j40
```

The project binary will be available at `build_release/art_modern`. Now we can test whether the program runs:

```shell
build_release/art_modern --help
build_release/art_modern --version # For version information
```

### Simulating WGS Data using _E. Coli_ Genome

Download _E. Coli_ reference genome from NCBI. Here we'll use K12 strand MG1655 substrand as an example.

```shell
mkdir -p tutorial_data
wget -4 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Now we can simulate WGS data using E. Coli reference genome. Let's satrt with single-end sequencing using HiSeq 2500 with 125bp read length and 10X coverage.

```shell
build_release/art_modern \
   --mode wgs \
   --lc se \
   --i-file tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna \
   --o-fastq tutorial_data/e_coli_wgs_se.fastq \
   --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
   --read_len 125 \
   --parallel 4 \
   --i-fcov 10
```

The generated FASTQ file will be at `tutorial_data/e_coli_wgs_se.fastq`.

We may also simulate paired-end data with following configuration:

```shell
build_release/art_modern \
   --mode wgs \
   --lc pe \
   --i-file tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna \
   --o-fastq tutorial_data/e_coli_wgs_pe.fastq \
   --qual_file_1 art/Illumina_profiles/HiSeq2500L125R1.txt \
   --qual_file_2 art/Illumina_profiles/HiSeq2500L125R2.txt \
   --read_len 125 \
   --parallel 4 \
   --i-fcov 10 \
   --pe_frag_dist_mean 300 \
   --pe_frag_dist_std_dev 50
```

Please note that we have additionally specified quality file for read 2 with the mean and standard diviation of fragment lengths.

### Simulating RNA-Seq Data using _C. Elegans_ Transcriptome

### Simulating Targeted Amplification Data

### Simulating WGS Data on Large Genomes

## Building

### Dependencies

- [CMake](https://cmake.org/), the building system used in our project. That further require:
  - A [CMake Generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html).
    - Under GNU/Linux, using [Ninja](https://ninja-build.org/) or [GNU Make](https://www.gnu.org/software/make) is preferred.
  - Under GNU/Linux, our project further requires [GNU Bash](https://www.gnu.org/software/bash).
- A working C and C++ compiler that supports C++ 17. The following compilers are supported:
  - [GCC](https://gcc.gnu.org/);
  - [Clang](https://clang.llvm.org/);
  - [Intel oneAPI DPCPP](https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/get-started-guide/2024-2/overview.html).
- [Boost C++ Library](https://www.boost.org/). This is an umbrella project of diverse small modules that can be used independently. Consult your system administrator for where to find them. The modules used in this project are namely:
  - **REQUIRED** [FileSystem](https://www.boost.org/doc/libs/1_85_0/libs/filesystem/);
  - **REQUIRED** [Regex](https://www.boost.org/doc/libs/1_85_0/libs/regex/);
  - **REQUIRED** [Program Options](https://www.boost.org/doc/libs/1_85_0/libs/program_options/);
  - **REQUIRED** [Thread](https://www.boost.org/doc/libs/1_85_0/libs/thread/);
  - **REQUIRED** [Log](https://www.boost.org/doc/libs/1_85_0/libs/log/);
  - **REQUIRED** [StackTrace](https://www.boost.org/doc/libs/1_85_0/doc/html/stacktrace.html);
  - **OPTIONAL** [Test](https://www.boost.org/doc/libs/1_85_0/libs/test/): For unit testing only. Can be absent for non-developers;
  - **OPTIONAL** [Timer](https://www.boost.org/doc/libs/1_85_0/libs/timer/): For displaying CPU and wall-clock time at the end of the program.
  - **OPTIONAL** `stacktrace_backtrace`: For a more developer-friendly stack trace. See [here](https://www.boost.org/doc/libs/1_85_0/doc/html/stacktrace/configuration_and_build.html) for details.
- A working [HTSLib](https://www.htslib.org/). You may either use the one bundled with the project or an external one that had already been installed inside your system.
  - To use bundled HTSLib sources, you need to have:
    - **REQUIRED** [zlib](https://www.zlib.net/);
    - **REQUIRED** [pthread](https://www.man7.org/linux/man-pages/man7/pthreads.7.html);
    - **OPTIONAL** [libbz2](http://www.bzip.org/): For CRAM compression, which is now not supported.
    - **OPTIONAL** [liblzma](https://tukaani.org/xz/): For CRAM compression, which is now not supported.
    - **OPTIONAL** [libdeflate](https://github.com/ebiggers/libdeflate): This library accelerates compressed BAM output.
    - See [official HTSLib documentation](https://github.com/samtools/samtools/blob/master/INSTALL) for more details. See also `USE_HTSLIB` CMake variable mentioned below.
  - To use external HTSLib, consult your system administrator. Those libraries usually named `libhts.so` with optional version suffixes.
- Optional libraries for accelerating random number generation. Including:
  - [Intel OneAPI Math Kernel Library (OneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html): Recommended for Intel or AMD CPUs.
  - [GNU Science Library (GSL)](https://www.gnu.org/software/gsl/): Used in original ART.
  - See `USE_RANDOM_GENERATOR` CMake variable mentioned below.

### CMake Variables

CMake variables control build behaviour. If you want a specific build (e.g., with accelerated random number generation, with or without debugging information), you should set them accordingly. They should be set when invoking `cmake`. For example,

```shell
cmake -DBUILD_SHARED_LIBS=ON
```

sets `BUILD_SHARED_LIBS` to `ON`.

- [`BUILD_SHARED_LIBS`](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html): Whether to build shared libraries.
  - **`ON` (DEFAULT): Will search for shared libraries and use dynamic linking.**
  - `OFF`: Will search for static libraries and use static linking.
- [`CMAKE_BUILD_TYPE`](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html): The CMake build type.
  - **`Debug` (DEFAULT): For developers with debugging needs.**
  - `Release`: Optimized executables/libraries without debug symbols.
  - `RelWithDebInfo`: Optimized executables/libraries with debug symbols.
- `CEU_CM_SHOULD_USE_NATIVE`: Whether to build the binaries using [`-mtune=native`](https://gcc.gnu.org/onlinedocs/gcc-14.1.0/gcc/x86-Options.html#index-march-16), if possible. This would result in faster executable with impaired portability (i.e., do not run on other machines).
  - **`OFF` (DEFAULT): Will not build native executables/libraries.**
  - `ON`: Will build native executables/libraries.
- `CEU_CM_SHOULD_ENABLE_TEST`: Whether test should be enabled.
  - **unset (DEFAULT): Set to `ON` if the CMake variable `CMAKE_BUILD_TYPE` is not `Release`, `OFF` otherwise.**
  - `OFF`: Will disable test.
  - `ON`: Will enable test.
- `USE_HTSLIB`: Use which HTSLib implementation.
  - **unset (DEFAULT): Will use bundled HTSLib.**
  - `hts`: Will use the HTSLib (`libhts.so`) found in the system.
  - Any other value `[val]`: Will use the HTSLib of other names (`lib[val].so`) found in the system.
- `USE_RANDOM_GENERATOR`: The random number generator used.
  - **`STL` (DEFAULT): Use STL random generators.**
  - `BOOST`: Use Boost random generators.
  - `GSL`: Use GSL random generators.
  - `ONEMKL`: Use Intel OneAPI MKL random generators. Highly recommended on Intel CPUs.
- `USE_THREAD_PARALLEL`: The thread-level parallelism strategy.
  - **`ASIO` (DEFAULT): Will use Boost ASIO for thread-based parallelism.**
  - `NOP`: Will not use thread-based parallelism.

## Usage

### Input Formats

Currently, we support input in FASTA and PBSim3 transcripts format.

**FOR FASTA FORMAT**: For read names, only characters before blank space are read.

A compatibility matrix is as follows:

| Parser \ Mode | `wgs`     | `trans`                     | `templ`                     |
|---------------|-----------|-----------------------------|-----------------------------|
| `memory`      | FASTA     | FASTA \| PBSim3 Transcripts | FASTA \| PBSim3 Transcripts |
| `htslib`      | FASTA     | **ERROR**                   | **ERROR**                   |
| `stream`      | **ERROR** | FASTA \| PBSim3 Transcripts | FASTA \| PBSim3 Transcripts |

### Library Construction Methods

### FASTA Parsers


## Acknowledgements

This simulator is based on the works of [Weichun Huang](mailto:whduke@gmail.com) _et al._, under [GNU GPL v3](https://www.gnu.org/licenses/) license. The software is originally distributed [here](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) with the following reference:

- W. Huang, L. Li, J. R. Myers, and G. T. Marth, _ART: a next-generation sequencing read simulator_, Bioinformatics (Oxford, England), vol. 28, no. 4, pp. 593–594, Feb. 2012, doi: [10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708).

The bundled HTSLib library used MIT License with the following reference:

- J. K. Bonfield et al., _HTSlib: C library for reading/writing high-throughput sequencing data_, GigaScience, vol. 10, no. 2, p. giab007, Jan. 2021, doi: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007).

## FAQ

### Supported CPUs?

Although this application should theoretically support both endianness, only little endian is tested. That is, if you're working on an Intel or AMD CPU, this application should work fine.

### How to split produced pair-end/mate-pair sequencing results to 2 FASTQ files?

This can be done through [`seqtk`](https://github.com/lh3/seqtk). For example, to split `tmp/test_small_pe/NC_001416.1.fq`:

```shell
# Read 1
seqtk seq tmp/test_small_pe/NC_001416.1.fq -1 > tmp/test_small_pe/NC_001416.1_1.fq
# Read 2
seqtk seq tmp/test_small_pe/NC_001416.1.fq -2 > tmp/test_small_pe/NC_001416.1_2.fq
```

You may also generate SAM/BAM files and extract PE FASTQ from them using `samtools`.

### Is there an interface for Python, Java, and more?

We may develop a C interface in the future after the APIs and design of the core library are settled.

### How to add adaptors \& primers to the reads?

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adapters to them, and use `templ` mode to introduce sequencing errors.

### I don't care about portability. How to make it wicked fast?

Easy. Set `USE_HTSLIB` to the latest HTSLib available on your system, `CMAKE_BUILD_TYPE` to `Release` or `RelWithDebInfo`, and `USE_RANDOM_GENERATOR` to `ONEMKL` on Intel/AMD machines. Please also make sure that your HTSLib has been linked with [libdeflate](https://github.com/ebiggers/libdeflate).

Also, when executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.
