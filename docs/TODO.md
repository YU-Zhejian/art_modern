# TODO

## IMPORTANT

- Some non-reproducible SEGFAULTs occured in `make testbuild` with pure-LLVM toolchain. CONFIG CMDLINE: 

  ```shell
  LD_LIBRARY_PATH=${HOME}/opt/boost-1.89.0-clang/lib/ \
      LD_RUN_PATH=${HOME}/opt/boost-1.89.0-clang/lib/ \
      PKG_CONFIG_PATH=${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/ \
      CMAKE_TOOLCHAIN_FILE=sh.d/toolchain/llvm-toolchain.cmake \
      /home/yuzj/miniconda3/envs/art_modern/bin/cmake \
      -G Ninja -Wdev -Wdeprecated --warn-uninitialized \
      -DCEU_CM_SHOULD_ENABLE_TEST=ON -DCMAKE_VERBOSE_MAKEFILE=ON \
      -DBoost_DIR=/home/yuzj/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/ \
      -DCMAKE_BUILD_TYPE=Debug \
      -DUSE_RANDOM_GENERATOR=STL \
      -DUSE_QUAL_GEN=WALKER \
      -DUSE_MALLOC=AUTO \
      -DUSE_THREAD_PARALLEL=ASIO \
      -DCMAKE_INSTALL_PREFIX=/tmp/art_modern-test-build-install-8faz6gpf \
      /home/yuzj/Documents/pbsim3_modern
  ```
  
  Failed output:

  ```text
  EXEC 68: /tmp/art_modern-test-build-install-8faz6gpf/bin/art_modern --builtin_qual_file HiSeq2500_125bp --i-file /home/yuzj/Documents/pbsim3_modern/data/raw_data/ce11.mRNA_head.fa --read_len 125 --i-batch_size 100 --mode template --lc pe --i-parser stream --i-fcov 10 --parallel 2 --ins_rate_1 0.1 --del_rate_1 0.1 --o-hl_sam /tmp/art_modern_test_small.d.TLjw30/test_small_pe_template_stream.hl.sam --o-fastq /tmp/art_modern_test_small.d.TLjw30/test_small_pe_template_stream.fq
  /home/yuzj/Documents/pbsim3_modern/sh.d/test-small.sh: line 109: 348299 Aborted                 (core dumped) env "ART_LOG_DIR=${OUT_DIR}/log_${EXEC_ORDER}.d" "${ART_CMD_ASSEMBLED[@]}" "$@" &>> "${OUT_DIR}"/am_exec_"${EXEC_ORDER}".log
  ```

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
- Pack MPI-enabled Debian packages, BioConda build, etc.
- Link MKL using pkgconfig using `mkl-sdl.pc` (Shipped with Intel) or `mkl-sdl-lp64.pc` (Shipped with Debian).

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams that allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worsen if POSIX AIO is used.

## I/O Formats

- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Support UCSC MAF output format?
- Add flags to disable/enable diverse BAM tags.

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
