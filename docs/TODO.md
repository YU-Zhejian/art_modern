# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- As BAM cannot accommodate contigs larger than 2GiB, should check this when initializing BAM writer. Suggest users to use unaligned BAM instead.
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.
- Add support for PE FASTA/FASTQ.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA, when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams, which allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worse if POSIX AIO is used.

## I/O Formats

- Support circular genome or RNA?
  - Search for new discussions in HTSLib, hts-specs, etc.
- Add `--i-nreads` to accurately specify the number of reads to simulate?

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
