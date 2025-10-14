# TODO

## IMPORTANT

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

- Supress all lintian issues.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams that allows reusing.
  - The current implementation passes too many small objects accross the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worsen if POSIX AIO is used.

- Support MPI-based parallelization. Basic ideas:
  - For `htslib` parser, just divide sequencing depth.
  - For `memory` parser, skip records based on MPI rank.
  - For `stream` parser, skip records based on MPI rank.
  - Revised dependencies section:
    - Optional MPI library for MPI-based parallelism. The following MPI implementations are supported:
      - [MPICH](https://www.mpich.org/).
      - [OpenMPI](https://www.open-mpi.org/).
      - [Intel MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html).

## I/O Formats

- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Support UCSC MAF output format?
- Add flags to disable/enable diverse BAM tags.

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
