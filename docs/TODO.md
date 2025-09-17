# TODO

## IMPORTANT

- ART stuck at one read. Why?
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.
- Hack Sphinx and Sphinx SVG-to-PDF converters to allow shields.io badges in the documentation like:

  ```markdown
  [![GitHub Release](https://img.shields.io/github/v/release/YU-Zhejian/art_modern)](https://github.com/YU-Zhejian/art_modern/)
  [![GitHub Downloads](https://img.shields.io/github/downloads/YU-Zhejian/art_modern/total.svg)](https://github.com/YU-Zhejian/art_modern/releases/)
  [![License](https://img.shields.io/badge/licence-GPL_3.0-blue.svg)](https://www.gnu.org/licenses/)
  
  [![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/art_modern/README.html)
  [![Conda Version](https://img.shields.io/conda/vn/bioconda/art_modern)](https://anaconda.org/bioconda/art_modern)
  [![Conda Downloads](https://img.shields.io/conda/dn/bioconda/art_modern)](https://anaconda.org/bioconda/art_modern)
  ```

## Packing

- Add copyright to all .cc/.hh files using one of the tools in <https://wiki.debian.org/CopyrightReviewTools>.
- Add manual pages.
- Generate RPM packages. See <https://rpm-packaging-guide.github.io/>.
- Is it a good idea to add `libpcg-cpp-dev` as compile dependency?
- Docker-ize the synthesis of DEBs, RPMs, Alpine Linux tarballs, etc.

## Performance

- The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams that allows reusing.
- Support MPI-based parallelization. Basic ideas:
  - For `htslib` parser, just divide sequencing depth.
  - For `memory` parser, skip records based on MPI rank.
  - For `stream` parser, skip records based on MPI rank.
  - Revised dependencies section:
    - Optional MPI library for MPI-based parallelism. The following MPI implementations are supported:
      - [MPICH](https://www.mpich.org/).
      - [OpenMPI](https://www.open-mpi.org/).
      - [Intel MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html).
  - Share the arguments between the main thread and the worker threads using pure MPI communication.

## I/O Formats

- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Support UCSC MAF output format?
- Support UCSC 2bit input format for fast on-disk random access of reference genome?
- Add flags to disable/enable diverse BAM tags.

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
