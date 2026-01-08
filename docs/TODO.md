# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- Perform large contig, many contig, and ultra deep tests.
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.
- `CEU_CM_SHOULD_ENABLE_TEST` seems not set by default. Correct it.
- Find concurrent queue and fmt using CMake.
- Supress unused args.
- Change year to 2026.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams that allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worsen if POSIX AIO is used.

## I/O Formats

- Support circular genome or RNA?
- Support simulating BGI/MGISEQ reads and new Illumina models: Started.
- Add `--i-nreads` to accurately specify the number of reads to simulate?

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
