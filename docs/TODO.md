# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- Perform large contig, many contig, and ultra deep tests.
- Add `art_modern-onenmpi` Conda package to Readme and News if approved.
- Support different read length for R1 and R2. Add testings. Update docs.
- Some non-reproducible SEGFAULTs occurred in `make testbuild` with pure-LLVM toolchain. CONFIG CMDLINE:

  ```shell
  LD_LIBRARY_PATH=${HOME}/opt/boost-1.89.0-clang/lib/ \
      LD_RUN_PATH=${HOME}/opt/boost-1.89.0-clang/lib/ \
      PKG_CONFIG_PATH=${HOME}/opt/fmt-12.0.0-clang/lib/pkgconfig/ \
      CMAKE_TOOLCHAIN_FILE=sh.d/toolchain/host-llvm/llvm-toolchain.cmake \
      /home/yuzj/miniconda3/envs/art_modern/bin/cmake \
      -G Ninja -Wdev -Wdeprecated --warn-uninitialized \
      -DCEU_CM_SHOULD_ENABLE_TEST=ON -DCMAKE_VERBOSE_MAKEFILE=ON \
      -DBoost_DIR=/home/yuzj/opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/ \
      -DCMAKE_BUILD_TYPE=Debug \
      -DUSE_RANDOM_GENERATOR=STL \
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

- Suppress all lintian issues.
- More testings required.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams that allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worsen if POSIX AIO is used.

## I/O Formats

- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Support UCSC MAF output format?
- Add flags to disable/enable diverse BAM tags.
- Support circular genome or RNA?
- Support simulating BGI/MGISEQ reads?
- Add `--i-nreads` to accurately specify the number of reads to simulate?

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.

## ERROR

```shell
/home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/bin/art_modern --i-file /home/yuzj/Documents/pbsim3_modern/data/raw_data/ce11.mRNA_head.pbsim3.transcript --i-type pbsim3_transcripts --i-batch_size 100 --mode template --lc pe --i-parser stream --parallel 4 --read_len_1 10 --read_len_2 150 --ins_rate_1 0.1 --del_rate_1 0.1 --o-hl_sam /tmp/art_modern_test_small.d.gzZZiN/test_small_pe_template_stream_pbsim3.hl.sam --o-fastq /tmp/art_modern_test_small.d.gzZZiN/test_small_pe_template_stream_pbsim3.fq
```

generates:

```text
[2025-11-11 22:33:51.729735] [T=0x000070a6e4e8d640] info: YuZJ Modified ART_Illumina (art_modern) v. 1.3.0 at <https://github.com/YU-Zhejian/art_modern/>
[2025-11-11 22:33:51.729817] [T=0x000070a6e4e8d640] info: Based on ART_Illumina: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)
[2025-11-11 22:33:51.729824] [T=0x000070a6e4e8d640] info: Originally written by: Weichun Huang <whduke@gmail.com>
[2025-11-11 22:33:51.729827] [T=0x000070a6e4e8d640] info: Modified by: YU Zhejian <yuzj25@seas.upenn.edu>
[2025-11-11 22:33:51.729829] [T=0x000070a6e4e8d640] info: Debugging functions enabled.
[2025-11-11 22:33:51.729907] [T=0x000070a6e4e8d640] info: Log file sink to '/tmp/art_modern_test_small.d.gzZZiN/log_22.d/nompi.log' added.
[2025-11-11 22:33:51.729945] [T=0x000070a6e4e8d640] info: MPI not found! Cross-node parallelism disabled.
[2025-11-11 22:33:51.730052] [T=0x000070a6e4e8d640] info: ARGS: /home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/bin/art_modern --i-file /home/yuzj/Documents/pbsim3_modern/data/raw_data/ce11.mRNA_head.pbsim3.transcript --i-type pbsim3_transcripts --i-batch_size 100 --mode template --lc pe --i-parser stream --parallel 4 --read_len_1 10 --read_len_2 150 --ins_rate_1 0.1 --del_rate_1 0.1 --o-hl_sam /tmp/art_modern_test_small.d.gzZZiN/test_small_pe_template_stream_pbsim3.hl.sam --o-fastq /tmp/art_modern_test_small.d.gzZZiN/test_small_pe_template_stream_pbsim3.fq
[2025-11-11 22:33:51.787876] [T=0x000070a6e4e8d640] info: QRange for R1: [3, 41].
[2025-11-11 22:33:51.788546] [T=0x000070a6e4e8d640] info: QRange for R2: [3, 41].
[2025-11-11 22:33:51.788554] [T=0x000070a6e4e8d640] info: Read quality profile loaded successfully.
[2025-11-11 22:33:51.788559] [T=0x000070a6e4e8d640] info: Read quality profile size for R1: 150
[2025-11-11 22:33:51.788565] [T=0x000070a6e4e8d640] info: Read quality profile size for R2: 150
[2025-11-11 22:33:51.789868] [T=0x000070a6e4e8d640] info: Argument parsing finished. Start generating...
[2025-11-11 22:33:51.789875] [T=0x000070a6e4e8d640] info: Boost::timer started.
[2025-11-11 22:33:51.794022] [T=0x000070a6e4e8d640] info: FASTQ LockFreeIO: Writer to '/tmp/art_modern_test_small.d.gzZZiN/test_small_pe_template_stream_pbsim3.fq' added.
[2025-11-11 22:33:51.834146] [T=0x000070a6e4e8d640] info: All writers added
[2025-11-11 22:33:51.840751] [T=0x000070a6e3d7f6c0] info: Starting simulation for job 1 with 100 contigs
[2025-11-11 22:33:51.845300] [T=0x000070a6e3d7f6c0] fatal: InMemoryFastaFetch::fetch: Requested range [4317, 4328) is out of bounds for sequence of length 4327.
[2025-11-11 22:33:51.845344] [T=0x000070a6e3d7f6c0] info: ABORT
[2025-11-11 22:33:51.851627] [T=0x000070a6e2d7d6c0] info: Starting simulation for job 2 with 100 contigs
[2025-11-11 22:33:51.854289] [T=0x000070a6e2d7d6c0] fatal: InMemoryFastaFetch::fetch: Requested range [2873, 2884) is out of bounds for sequence of length 2883.
[2025-11-11 22:33:51.854368] [T=0x000070a6e2d7d6c0] info: ABORT
[2025-11-11 22:33:51.862778] [T=0x000070a6e357e6c0] info: Starting simulation for job 3 with 100 contigs
[2025-11-11 22:33:51.877004] [T=0x000070a6e357e6c0] fatal: InMemoryFastaFetch::fetch: Requested range [955, 966) is out of bounds for sequence of length 965.
[2025-11-11 22:33:51.877039] [T=0x000070a6e357e6c0] info: ABORT
[2025-11-11 22:33:51.882169] [T=0x000070a6e257c6c0] info: Starting simulation for job 4 with 100 contigs
[2025-11-11 22:33:51.906066] [T=0x000070a6e257c6c0] fatal: InMemoryFastaFetch::fetch: Requested range [583, 595) is out of bounds for sequence of length 593.
[2025-11-11 22:33:51.906103] [T=0x000070a6e257c6c0] info: ABORT
[2025-11-11 22:33:51.948186] [T=0x000070a6e4e8d640] info: All jobs submitted. Waiting for job pool to stop...
[2025-11-11 22:33:51.845361] [T=0x000070a6e3d7f6c0] info: Stacktrace:
 0# labw::art_modern::abort_mpi(int) at /home/yuzj/Documents/pbsim3_modern/src/libam_support/utils/mpi_utils.cc:90
 1# labw::art_modern::InMemoryFastaFetch::fetch[abi:cxx11](unsigned long, long, long) at /home/yuzj/Documents/pbsim3_modern/src/libam_support/ref/fetch/InMemoryFastaFetch.cc:98
 2# labw::art_modern::ArtContig::generate_read_pe(bool, labw::art_modern::ArtRead&, labw::art_modern::ArtRead&) at /home/yuzj/Documents/pbsim3_modern/src/art/lib/ArtContig.cc:100
 3# labw::art_modern::ArtJobExecutor::generate_pe_(labw::art_modern::ArtContig&, bool, long, labw::art_modern::Rprob&) at /home/yuzj/Documents/pbsim3_modern/src/art/lib/ArtJobExecutor.cc:129
 4# labw::art_modern::ArtJobExecutor::generate_(long, bool, labw::art_modern::ArtContig&, long&, labw::art_modern::Rprob&) at /home/yuzj/Documents/pbsim3_modern/src/art/lib/ArtJobExecutor.cc:102
 5# labw::art_modern::ArtJobExecutor::operator()() at /home/yuzj/Documents/pbsim3_modern/src/art/lib/ArtJobExecutor.cc:230
 6# boost::asio::detail::executor_op<boost::asio::detail::binder0<labw::art_modern::JobPool::add(std::shared_ptr<labw::art_modern::JobExecutor> const&)::{lambda()#1}>, std::allocator<void>, boost::asio::detail::scheduler_operation>::do_complete(void*, boost::asio::detail::scheduler_operation*, boost::system::error_code const&, unsigned long) at /usr/include/boost/asio/detail/executor_op.hpp:74
 7# boost::asio::detail::scheduler::do_run_one(boost::asio::detail::conditionally_enabled_mutex::scoped_lock&, boost::asio::detail::scheduler_thread_info&, boost::system::error_code const&) in /home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/lib/art_modern/lib/libam_support_lib.so
 8# boost::asio::detail::scheduler::run(boost::system::error_code&) in /home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/lib/art_modern/lib/libam_support_lib.so
 9# boost::asio::thread_pool::thread_function::operator()() in /home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/lib/art_modern/lib/libam_support_lib.so
10# boost::asio::detail::posix_thread::func<boost::asio::thread_pool::thread_function>::run() in /home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/lib/art_modern/lib/libam_support_lib.so
11# boost_asio_detail_posix_thread_function in /home/yuzj/Documents/pbsim3_modern/opt/build_debug_install/lib/art_modern/lib/libam_support_lib.so
12# start_thread at ./nptl/pthread_create.c:447
13# clone3 at ../sysdeps/unix/sysv/linux/x86_64/clone3.S:80
```