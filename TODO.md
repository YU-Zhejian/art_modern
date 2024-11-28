# TODO

- Make it even faster.
- Support MPI-based parallelization. Basic ideas:
  - For `htslib` parser, just divide sequencing depth.
  - For `memory` parser, skip records based on MPI rank.
  - For `stream` parser, skip records based on MPI rank.
  - Revised dependencies section:
    - Optional MPI library for MPI-based parallelism.
      - MPI implementations (library and compiler). The following MPI implementations are supported:
        - [MPICH](https://www.mpich.org/).
        - [OpenMPI](https://www.open-mpi.org/).
        - [Intel MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html).
        - [MS-MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi) (For working under MSYS2). See also: [MSYS2 Package Repository](https://packages.msys2.org/packages/mingw-w64-x86_64-msmpi)
      - Google Protocol Buffers (Protobuf) for serialization/deserialization of MPI.
- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Support UCSC MAF output format?
- Working with >65535 CIGAR operations? See <https://github.com/lh3/minimap2?tab=readme-ov-file#working-with-65535-cigar-operations>
- Support UCSC 2bit input format?

## Random Number Generation Functions

The current random number generation function in each library is MT19937, which may not be the best choice for performance-critical applications. However, it is the most widely used, well-known, and is implemented in all random number generator libraries (namely, Boost, GSL, STL, and Intel OneAPI MKL). We may further introduce Taus RNGs in the future.

