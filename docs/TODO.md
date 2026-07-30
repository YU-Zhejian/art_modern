# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.
- Add support for PE FASTA/FASTQ.

## Possible Simplifications

- Drop `USE_MALLOC` from CMake entirely. jemalloc, mimalloc and tcmalloc can all be `LD_PRELOAD`ed at
  run time, so selecting one at link time buys little for a shared build. Removing it would delete
  `src/cmake.d/malloc.cmake` (146 lines), the `USE_MALLOC` axis of the `test-build.py` matrix and its
  `lib{je,mi,tc}malloc*.so` link assertions, the `USE_MALLOC` section of `docs/Install.md`, and the
  `-DUSE_MALLOC=NOP` in `Makefile` and `packing/pack-alpine-exe.sh`.
  - Caveats to weigh first: `LD_PRELOAD` does not work for the fully static Alpine build, though that
    one already uses `NOP`; and under MPI the variable has to reach every rank, which depends on the
    launcher forwarding the environment.

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
