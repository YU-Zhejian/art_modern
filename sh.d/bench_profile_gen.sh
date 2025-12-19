#!/usr/bin/env bash

mkdir -p opt/bench_profile_gen
opt/build_release_install/bin/art_modern \
    --i-file data/raw_data/ce11_chr1.fa --i-fcov 400 \
    --o-sam opt/bench_profile_gen/out_pe.bam \
    --o-sam-write_bam \
    --o-fastq opt/bench_profile_gen/out_pe.fq \
    --pe_frag_dist_mean 200 \
    --pe_frag_dist_std_dev 50 \
    --lc pe
# Order: Old, new
#[2025-12-18 21:29:52.124864] [T=0x00007099780ee640] info: Total reads processed: 40.19M speed: 855.42K reads/s
#[2025-12-18 21:29:52.124891] [T=0x00007099780ee640] info: Total bases processed: 6.03G speed: 128.31M bases/s
#[2025-12-18 21:29:52.124900] [T=0x00007099780ee640] info: Time elapsed:  46.990000s wall, 767.130000s user + 53.950000s system = 821.080000s CPU (1747.4%)
#[2025-12-18 21:26:23.579286] [T=0x0000777414efd640] info: Total reads processed: 40.19M speed: 4.96M reads/s
#[2025-12-18 21:26:23.579304] [T=0x0000777414efd640] info: Total bases processed: 6.03G speed: 743.50M bases/s
#[2025-12-18 21:26:23.579311] [T=0x0000777414efd640] info: Time elapsed:  8.110000s wall, 70.460000s user + 2.800000s system = 73.260000s CPU (903.3%)
opt/build_release_install/bin/art_profile_builder \
    --parallel 20 \
    --o-file1 /dev/null \
    --o-file2 /dev/null \
    --is_pe \
    --i-file opt/bench_profile_gen/out_pe.bam \
    --read_len 150
#[2025-12-18 21:31:18.236044] [T=0x000077feb670f640] info: Total reads processed: 40.19M speed: 1.44M reads/s
#[2025-12-18 21:31:18.236073] [T=0x000077feb670f640] info: Total bases processed: 6.03G speed: 215.46M bases/s
#[2025-12-18 21:31:18.236080] [T=0x000077feb670f640] info: Time elapsed:  27.980000s wall, 386.830000s user + 79.380000s system = 466.210000s CPU (1666.2%)
#[2025-12-18 21:27:10.544972] [T=0x00007da25a4f7640] info: Total reads processed: 40.19M speed: 1.61M reads/s
#[2025-12-18 21:27:10.544997] [T=0x00007da25a4f7640] info: Total bases processed: 6.03G speed: 241.60M bases/s
#[2025-12-18 21:27:10.545008] [T=0x00007da25a4f7640] info: Time elapsed:  24.950000s wall, 60.710000s user + 8.090000s system = 68.800000s CPU (275.8%)
opt/build_release_install/bin/art_profile_builder \
    --parallel 20 \
    --o-file1 /dev/null \
    --o-file2 /dev/null \
    --is_pe \
    --i-file opt/bench_profile_gen/out_pe.fq \
    --read_len 150
