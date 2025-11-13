# Readme

## Many Contigs

Many contig stress test. Simulators that use integer to represent contig ID will fail.

Here, a generator will generate 5G contigs (`contig000000000` to `contig140000000`), each 1K bases long and 5T bases in total, to `/dev/stdout`.

## Large Contigs

Large contigs stress test. Simulators that use integer to represent contig offset will fail.

Here, a generator will generate 4 contigs (`contig000000000` to `contig000000004`), each 5G bases long and 20G bases in total, to `/dev/stdout`. Under a 20X sequence coverage and SE 125, the simulator should generate 3.2G reads.

Note that this test organism have smaller genome than _Lepidosiren paradoxa_ (South American lungfish) [ASM4058144v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040581445.1/).
