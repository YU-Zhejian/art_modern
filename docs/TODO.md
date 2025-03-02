# TODO

## IMPORTANT

- ART stuck at one read. Why?
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.

## Performance

- The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. Consider using the method implemented in `pigz`.
- The current implementation tightly couples Moody Camel queue with output writers, which is not wise.
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
- Revise support over other random number generation functions.
- Builtin profiles take too much space on the executable. May consider:
  - Use a CMake option that disables embedding of builtin profiles.

## I/O Formats

- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Support UCSC MAF output format?
- Support UCSC 2bit input format for fast on-disk random access of reference genome?
- Add flags to disable/enable diverse BAM tags.

## Known Incompatibilities

- GCC 13.3.0 on Haiku OS hrev58590 may generate a kernel panic that jam the entire system.
- GCC would fail on Debian GNU/Hurd.
