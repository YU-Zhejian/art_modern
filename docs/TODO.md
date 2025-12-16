# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- Perform large contig, many contig, and ultra deep tests.
- The `art_modern` would still stuck on random situations.

  Under Intel compiler configuration, with args:

  ```text
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_BUILD_TYPE=Debug -DUSE_RANDOM_GENERATOR=PCG -DUSE_MALLOC=AUTO -DUSE_THREAD_PARALLEL=ASIO -DAM_NO_Q_REVERSE=ON -DAM_NO_Q_REVERSE=ON
  ```

  ```text
  [2025-12-12 18:31:02.372779] [T=0x00007d21a70fe1c0] info: YuZJ Modified ART_Illumina (art_modern) v. 1.3.1 at <https://github.com/YU-Zhejian/art_modern/>
  [2025-12-12 18:31:02.372844] [T=0x00007d21a70fe1c0] info: Based on ART_Illumina: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)
  [2025-12-12 18:31:02.372850] [T=0x00007d21a70fe1c0] info: Originally written by: Weichun Huang <whduke@gmail.com>
  [2025-12-12 18:31:02.372852] [T=0x00007d21a70fe1c0] info: Modified by: YU Zhejian <yuzj25@seas.upenn.edu>
  [2025-12-12 18:31:02.372854] [T=0x00007d21a70fe1c0] info: Debugging functions enabled.
  [2025-12-12 18:31:02.372861] [T=0x00007d21a70fe1c0] info: ART_LOG_DIR defined; Using '<PROJDIR>test-build-6476.log.d/test-small-out-11/log_102.d' as log directory.
  [2025-12-12 18:31:02.372928] [T=0x00007d21a70fe1c0] info: Log file sink to '<PROJDIR>test-build-6476.log.d/test-small-out-11/log_102.d/nompi.log' added.
  [2025-12-12 18:31:02.372959] [T=0x00007d21a70fe1c0] info: MPI not found! Cross-node parallelism disabled.
  [2025-12-12 18:31:02.373036] [T=0x00007d21a70fe1c0] info: ARGS: <PROJDIR>test-build-6476.log.d/art_modern-test-install-11/bin/art_modern --i-file <PROJDIR>/data/raw_data/ce11.mRNA_head.pbsim3.transcript --i-type pbsim3_transcripts --i-batch_size 100 --mode template --lc pe --i-parser stream --parallel 2 --read_len_1 10 --read_len_2 150 --ins_rate_1 0.1 --del_rate_1 0.1 --o-hl_sam <PROJDIR>test-build-6476.log.d/test-small-out-11/test_small_pe_template_stream_pbsim3.hl.sam --o-fastq <PROJDIR>test-build-6476.log.d/test-small-out-11/test_small_pe_template_stream_pbsim3.fq
  [2025-12-12 18:31:02.406309] [T=0x00007d21a70fe1c0] info: QRange for R1: [3, 41].
  [2025-12-12 18:31:02.406790] [T=0x00007d21a70fe1c0] info: QRange for R2: [3, 41].
  [2025-12-12 18:31:02.406797] [T=0x00007d21a70fe1c0] info: Read quality profile loaded successfully.
  [2025-12-12 18:31:02.406801] [T=0x00007d21a70fe1c0] info: Read quality profile size for R1: 150
  [2025-12-12 18:31:02.406806] [T=0x00007d21a70fe1c0] info: Read quality profile size for R2: 150
  [2025-12-12 18:31:02.407372] [T=0x00007d21a70fe1c0] info: Argument parsing finished. Start generating...
  [2025-12-12 18:31:02.407376] [T=0x00007d21a70fe1c0] info: Boost::timer started.
  [2025-12-12 18:31:02.410482] [T=0x00007d21a70fe1c0] info: FASTQ LockFreeIO: Writer to '<PROJDIR>test-build-6476.log.d/test-small-out-11/test_small_pe_template_stream_pbsim3.fq' added.
  [2025-12-12 18:31:02.432878] [T=0x00007d21a70fe1c0] info: All writers added
  [2025-12-12 18:31:02.433526] [T=0x00007d21a3f7f6c0] info: Starting simulation for job 1 with 100 contigs
  [2025-12-12 18:31:02.433872] [T=0x00007d21a377e6c0] info: Starting simulation for job 2 with 100 contigs
  [2025-12-12 18:31:02.437800] [T=0x00007d21a70fe1c0] info: All jobs submitted. Waiting for job pool to stop...
  [2025-12-12 18:31:02.454294] [T=0x00007d21a377e6c0] info: Finished simulation for job 2 with 510.00 reads (mean depth=510) generated.
  [2025-12-12 18:31:02.454307] [T=0x00007d21a3f7f6c0] info: Finished simulation for job 1 with 538.00 reads (mean depth=538) generated.
  [2025-12-12 18:31:03.433619] [T=0x00007d21a3f7f6c0] info: Starting simulation for job 3 with 100 contigs
  [2025-12-12 18:31:03.433916] [T=0x00007d21a377e6c0] info: Starting simulation for job 4 with 100 contigs
  [2025-12-12 18:31:03.447323] [T=0x00007d21a377e6c0] info: Finished simulation for job 4 with 434.00 reads (mean depth=434) generated.
  [2025-12-12 18:31:03.447766] [T=0x00007d21a3f7f6c0] info: Finished simulation for job 3 with 484.00 reads (mean depth=484) generated.
  [2025-12-12 18:31:04.433784] [T=0x00007d21a3f7f6c0] info: Starting simulation for job 5 with 100 contigs
  [2025-12-12 18:31:04.434059] [T=0x00007d21a377e6c0] info: Starting simulation for job 6 with 100 contigs
  [2025-12-12 18:31:04.449879] [T=0x00007d21a3f7f6c0] info: Finished simulation for job 5 with 498.00 reads (mean depth=498) generated.
  [2025-12-12 18:31:04.450382] [T=0x00007d21a377e6c0] info: Finished simulation for job 6 with 524.00 reads (mean depth=524) generated.
  [2025-12-12 18:31:05.434042] [T=0x00007d21a3f7f6c0] info: Starting simulation for job 7 with 100 contigs
  [2025-12-12 18:31:05.434287] [T=0x00007d21a377e6c0] info: Starting simulation for job 8 with 100 contigs
  [2025-12-12 18:31:05.486728] [T=0x00007d21a3f7f6c0] info: Finished simulation for job 7 with 518.00 reads (mean depth=518) generated.
  [2025-12-12 18:31:06.434172] [T=0x00007d21a3f7f6c0] info: Starting simulation for job 9 with 100 contigs
  [2025-12-12 18:31:06.456935] [T=0x00007d21a3f7f6c0] info: Finished simulation for job 9 with 490.00 reads (mean depth=490) generated.
  [2025-12-12 18:31:07.434325] [T=0x00007d21a3f7f6c0] info: Starting simulation for job 10 with 100 contigs
  [2025-12-12 18:31:07.453710] [T=0x00007d21a3f7f6c0] info: Finished simulation for job 10 with 472.00 reads (mean depth=472) generated.
  [2025-12-12 18:31:10.434675] [T=0x00007d2197fff6c0] info: AJEReporter: Job 8:nompi | ON: 'NM_062251' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 120.00
  [2025-12-12 18:31:12.408529] [T=0x00007d21a277c6c0] info: JobPoolReporter: 0 JobExecutors running
  [2025-12-12 18:31:15.436637] [T=0x00007d2197fff6c0] info: AJEReporter: Job 8:nompi | ON: 'NM_062251' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 120.00
  [2025-12-12 18:31:20.436976] [T=0x00007d2197fff6c0] info: AJEReporter: Job 8:nompi | ON: 'NM_062251' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 120.00
  [2025-12-12 18:31:22.410541] [T=0x00007d21a277c6c0] info: JobPoolReporter: 0 JobExecutors running
  [2025-12-12 18:31:25.437413] [T=0x00007d2197fff6c0] info: AJEReporter: Job 8:nompi | ON: 'NM_062251' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 120.00
  [2025-12-12 18:31:30.437750] [T=0x00007d2197fff6c0] info: AJEReporter: Job 8:nompi | ON: 'NM_062251' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 120.00
  [...]
  killed
  ```

  Same compiler, with args:

  ```text
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_BUILD_TYPE=Release -DUSE_RANDOM_GENERATOR=ONEMKL -DUSE_MALLOC=AUTO -DUSE_THREAD_PARALLEL=ASIO -DWITH_MPI=ON
  ```

  ```text
  [2025-12-12 20:16:29.804260] [T=0x000074a3ef62f200] info: YuZJ Modified ART_Illumina (art_modern) v. 1.3.1 at <https://github.com/YU-Zhejian/art_modern/>
  [2025-12-12 20:16:29.804401] [T=0x000074a3ef62f200] info: Based on ART_Illumina: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)
  [2025-12-12 20:16:29.804411] [T=0x000074a3ef62f200] info: Originally written by: Weichun Huang <whduke@gmail.com>
  [2025-12-12 20:16:29.804415] [T=0x000074a3ef62f200] info: Modified by: YU Zhejian <yuzj25@seas.upenn.edu>
  [2025-12-12 20:16:29.804426] [T=0x000074a3ef62f200] info: ART_LOG_DIR defined; Using '<PROJDIR>test-build-393427.log.d/test-small-out-33/log_102.d' as log directory.
  [2025-12-12 20:16:29.804534] [T=0x000074a3ef62f200] info: Log file sink to '<PROJDIR>test-build-393427.log.d/test-small-out-33/log_102.d/0.log' added.
  [2025-12-12 20:16:29.804574] [T=0x000074a3ef62f200] info: MPI found! Cross-node parallelism enabled.
  [2025-12-12 20:16:29.804582] [T=0x000074a3ef62f200] info: MPI main process started.
  [2025-12-12 20:16:29.804687] [T=0x000074a3ef62f200] info: ARGS: <PROJDIR>test-build-393427.log.d/art_modern-test-install-33/bin/art_modern-mpi --i-file <PROJDIR>/data/raw_data/ce11.mRNA_head.pbsim3.transcript --i-type pbsim3_transcripts --i-batch_size 100 --mode template --lc pe --i-parser stream --parallel 2 --read_len_1 10 --read_len_2 150 --ins_rate_1 0.1 --del_rate_1 0.1 --o-hl_sam <PROJDIR>test-build-393427.log.d/test-small-out-33/test_small_pe_template_stream_pbsim3.hl.sam --o-fastq <PROJDIR>test-build-393427.log.d/test-small-out-33/test_small_pe_template_stream_pbsim3.fq
  [2025-12-12 20:16:29.847016] [T=0x000074a3ef62f200] info: QRange for R1: [3, 41].
  [2025-12-12 20:16:29.847776] [T=0x000074a3ef62f200] info: QRange for R2: [3, 41].
  [2025-12-12 20:16:29.847797] [T=0x000074a3ef62f200] info: Read quality profile loaded successfully.
  [2025-12-12 20:16:29.847809] [T=0x000074a3ef62f200] info: Read quality profile size for R1: 150
  [2025-12-12 20:16:29.847823] [T=0x000074a3ef62f200] info: Read quality profile size for R2: 150
  [2025-12-12 20:16:29.848699] [T=0x000074a3ef62f200] info: Argument parsing finished. Start generating...
  [2025-12-12 20:16:29.848715] [T=0x000074a3ef62f200] info: Boost::timer started.
  [2025-12-12 20:16:29.854530] [T=0x000074a3ef62f200] info: FASTQ LockFreeIO: Writer to '<PROJDIR>test-build-393427.log.d/test-small-out-33/test_small_pe_template_stream_pbsim3.0.fq' added.
  [2025-12-12 20:16:29.891454] [T=0x000074a3ef62f200] info: All writers added
  [2025-12-12 20:16:29.897571] [T=0x000074a3ef62f200] info: All jobs submitted. Waiting for job pool to stop...
  [2025-12-12 20:16:29.913098] [T=0x000074a34a1fe6c0] info: Starting simulation for job 1 with 100 contigs
  [2025-12-12 20:16:29.913117] [T=0x000074a34a9ff6c0] info: Starting simulation for job 2 with 100 contigs
  [2025-12-12 20:16:29.920558] [T=0x000074a34a9ff6c0] info: Finished simulation for job 2 with 482.00 reads (mean depth=482) generated.
  [2025-12-12 20:16:30.893879] [T=0x000074a34a9ff6c0] info: Starting simulation for job 3 with 50 contigs
  [2025-12-12 20:16:30.897473] [T=0x000074a34a9ff6c0] info: Finished simulation for job 3 with 258.00 reads (mean depth=258) generated.
  [2025-12-12 20:16:34.892953] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:16:39.850709] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:16:39.896631] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:16:44.896960] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:16:49.851334] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:16:49.897419] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:16:54.900496] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:16:59.851979] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:16:59.900854] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:04.901280] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:09.853006] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:17:09.905433] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:14.906604] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:19.853647] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:17:19.909706] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:24.910134] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:29.854551] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:17:29.910491] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:34.910974] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [2025-12-12 20:17:39.856476] [T=0x000074a3491fc6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-12 20:17:39.911347] [T=0x000074a33f5ff6c0] info: AJEReporter: Job 1:0 | ON: 'NM_001027784' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 122.00
  [...]
  killed
  ```

  Clang compiler, with args:

  ```text
  -DBoost_DIR=<HOMEDIR>opt/boost-1.89.0-clang/lib/cmake/Boost-1.89.0/ -DCMAKE_BUILD_TYPE=Debug -DUSE_RANDOM_GENERATOR=PCG -DUSE_MALLOC=JEMALLOC -DUSE_THREAD_PARALLEL=ASIO
  ```

  ```text
  [2025-12-15 15:43:42.380385] [T=0x00007cb002ef61c0] info: YuZJ Modified ART_Illumina (art_modern) v. 1.3.1 at <https://github.com/YU-Zhejian/art_modern/>
  [2025-12-15 15:43:42.380514] [T=0x00007cb002ef61c0] info: Based on ART_Illumina: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)
  [2025-12-15 15:43:42.380521] [T=0x00007cb002ef61c0] info: Originally written by: Weichun Huang <whduke@gmail.com>
  [2025-12-15 15:43:42.380523] [T=0x00007cb002ef61c0] info: Modified by: YU Zhejian <yuzj25@seas.upenn.edu>
  [2025-12-15 15:43:42.380525] [T=0x00007cb002ef61c0] info: Debugging functions enabled.
  [2025-12-15 15:43:42.380532] [T=0x00007cb002ef61c0] info: ART_LOG_DIR defined; Using '<PROJDIR>test-build-2099768.log.d/test-small-out-4/log_102.d' as log directory.
  [2025-12-15 15:43:42.380637] [T=0x00007cb002ef61c0] info: Log file sink to '<PROJDIR>test-build-2099768.log.d/test-small-out-4/log_102.d/nompi.log' added.
  [2025-12-15 15:43:42.380664] [T=0x00007cb002ef61c0] info: MPI not found! Cross-node parallelism disabled.
  [2025-12-15 15:43:42.380743] [T=0x00007cb002ef61c0] info: ARGS: <PROJDIR>test-build-2099768.log.d/art_modern-test-install-4/bin/art_modern --i-file <PROJDIR>/data/raw_data/ce11.mRNA_head.pbsim3.transcript --i-type pbsim3_transcripts --i-batch_size 100 --mode template --lc pe --i-parser stream --parallel 2 --read_len_1 10 --read_len_2 150 --ins_rate_1 0.1 --del_rate_1 0.1 --o-hl_sam <PROJDIR>test-build-2099768.log.d/test-small-out-4/test_small_pe_template_stream_pbsim3.hl.sam --o-fastq <PROJDIR>test-build-2099768.log.d/test-small-out-4/test_small_pe_template_stream_pbsim3.fq
  [2025-12-15 15:43:42.464132] [T=0x00007cb002ef61c0] info: QRange for R1: [3, 41].
  [2025-12-15 15:43:42.465730] [T=0x00007cb002ef61c0] info: QRange for R2: [3, 41].
  [2025-12-15 15:43:42.465760] [T=0x00007cb002ef61c0] info: Read quality profile loaded successfully.
  [2025-12-15 15:43:42.465774] [T=0x00007cb002ef61c0] info: Read quality profile size for R1: 150
  [2025-12-15 15:43:42.465789] [T=0x00007cb002ef61c0] info: Read quality profile size for R2: 150
  [2025-12-15 15:43:42.466788] [T=0x00007cb002ef61c0] info: Argument parsing finished. Start generating...
  [2025-12-15 15:43:42.466812] [T=0x00007cb002ef61c0] info: Boost::timer started.
  [2025-12-15 15:43:42.472506] [T=0x00007cb002ef61c0] info: FASTQ LockFreeIO: Writer to '<PROJDIR>test-build-2099768.log.d/test-small-out-4/test_small_pe_template_stream_pbsim3.fq' added.
  [2025-12-15 15:43:42.507930] [T=0x00007cb002ef61c0] info: All writers added
  [2025-12-15 15:43:42.513700] [T=0x00007cafff37f6c0] info: Starting simulation for job 1 with 100 contigs
  [2025-12-15 15:43:42.515692] [T=0x00007caffeb7e6c0] info: Starting simulation for job 2 with 100 contigs
  [2025-12-15 15:43:42.527708] [T=0x00007cb002ef61c0] info: All jobs submitted. Waiting for job pool to stop...
  [2025-12-15 15:43:42.540258] [T=0x00007cafff37f6c0] info: Finished simulation for job 1 with 538.00 reads (mean depth=538) generated.
  [2025-12-15 15:43:42.540484] [T=0x00007cafff37f6c0] info: Starting simulation for job 3 with 100 contigs
  [2025-12-15 15:43:42.594480] [T=0x00007cafff37f6c0] info: Finished simulation for job 3 with 484.00 reads (mean depth=484) generated.
  [2025-12-15 15:43:42.594652] [T=0x00007cafff37f6c0] info: Starting simulation for job 4 with 100 contigs
  [2025-12-15 15:43:42.616864] [T=0x00007cafff37f6c0] info: Finished simulation for job 4 with 434.00 reads (mean depth=434) generated.
  [2025-12-15 15:43:42.616989] [T=0x00007cafff37f6c0] info: Starting simulation for job 5 with 100 contigs
  [2025-12-15 15:43:42.637322] [T=0x00007cafff37f6c0] info: Finished simulation for job 5 with 498.00 reads (mean depth=498) generated.
  [2025-12-15 15:43:42.637438] [T=0x00007cafff37f6c0] info: Starting simulation for job 6 with 100 contigs
  [2025-12-15 15:43:42.656891] [T=0x00007cafff37f6c0] info: Finished simulation for job 6 with 524.00 reads (mean depth=524) generated.
  [2025-12-15 15:43:42.657015] [T=0x00007cafff37f6c0] info: Starting simulation for job 7 with 100 contigs
  [2025-12-15 15:43:42.677934] [T=0x00007cafff37f6c0] info: Finished simulation for job 7 with 518.00 reads (mean depth=518) generated.
  [2025-12-15 15:43:42.678027] [T=0x00007cafff37f6c0] info: Starting simulation for job 8 with 100 contigs
  [2025-12-15 15:43:42.698325] [T=0x00007cafff37f6c0] info: Finished simulation for job 8 with 500.00 reads (mean depth=500) generated.
  [2025-12-15 15:43:42.698415] [T=0x00007cafff37f6c0] info: Starting simulation for job 9 with 100 contigs
  [2025-12-15 15:43:42.718479] [T=0x00007cafff37f6c0] info: Finished simulation for job 9 with 490.00 reads (mean depth=490) generated.
  [2025-12-15 15:43:42.718551] [T=0x00007cafff37f6c0] info: Starting simulation for job 10 with 100 contigs
  [2025-12-15 15:43:42.738613] [T=0x00007cafff37f6c0] info: Finished simulation for job 10 with 472.00 reads (mean depth=472) generated.
  [2025-12-15 15:43:47.516674] [T=0x00007caff1dff6c0] info: AJEReporter: Job 2:nompi | ON: 'NM_071427' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 92.00
  [2025-12-15 15:43:52.467528] [T=0x00007caffdb7c6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-15 15:43:52.516832] [T=0x00007caff1dff6c0] info: AJEReporter: Job 2:nompi | ON: 'NM_071427' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 92.00
  [2025-12-15 15:43:57.517005] [T=0x00007caff1dff6c0] info: AJEReporter: Job 2:nompi | ON: 'NM_071427' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 92.00
  [2025-12-15 15:44:02.467680] [T=0x00007caffdb7c6c0] info: JobPoolReporter: 1 JobExecutors running
  [2025-12-15 15:44:02.517163] [T=0x00007caff1dff6c0] info: AJEReporter: Job 2:nompi | ON: 'NM_071427' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 92.00
  [2025-12-15 15:44:07.525510] [T=0x00007caff1dff6c0] info: AJEReporter: Job 2:nompi | ON: 'NM_071427' | SUCCESS: current=2.00, left=2.00 | FAIL: current=0.00, max=5.00 | TOTAL: 92.00
  [2025-12-15 15:44:12.470488] [T=0x00007caffdb7c6c0] info: JobPoolReporter: 1 JobExecutors running
  [...]
  killed
  ```

- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.

## Packing

- Suppress all lintian issues.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams that allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worsen if POSIX AIO is used.

## I/O Formats

- Support [Illumina Complete Long Read](https://www.illumina.com/products/by-brand/complete-long-reads-portfolio.html)?
- Add flags to disable/enable diverse BAM tags.
- Support circular genome or RNA?
- Support simulating BGI/MGISEQ reads?
- Add `--i-nreads` to accurately specify the number of reads to simulate?

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.

## Simulate scRNA-Seq Data

1. Add gene body coverage bias introduced by scRNA-Seq platforms.
2. Add cell barcode and UMI bias.
3. Add gene strand bias.
4. The generated cell barcode and UMIs can be used as input to another ART simulation to add sequencing errors.
