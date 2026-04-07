# Readme

## Limitations of BAM format

1. BAM uses `uint32t` to record number of contigs and uses `int32_t` for contig ID, whose actual range should be -1 to $2^31-1$ inclusive.
2. BAI cannot handle contigs larger than 512Mbp. CSI may handle 2Gbp contigs.
3. BAM `POS` field is `int32_t`, whose actual range should be -1 to $2^31-1$ inclusive.

See [this](https://github.com/samtools/hts-specs/issues/240), [this](https://github.com/samtools/hts-specs/issues/40), and [this](https://github.com/samtools/hts-specs/issues/735) discussion in HTSLib-related repos.

## Many Contigs

Many contig stress test. Simulators that use `uint32_t` to represent contig ID will fail.

Here, a generator will generate 5G contigs (`contig000000000` to `contig140000000`), each 1K bases long and 5T bases in total, to `/dev/stdout`.

## Large Contigs

Large contigs stress test. Simulators that use `uint32_t` to represent contig offset will fail.

Here, a generator will generate 4 contigs (`contig000000000` to `contig000000004`), each 5G bases long and 20G bases in total, to `/dev/stdout`. Under a 20X sequence coverage and SE 125, the simulator should generate 3.2G reads.

Note that this test organism have smaller genome than _Lepidosiren paradoxa_ (South American lungfish) [ASM4058144v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040581445.1/).

The test would also check whether the produced FASTQ read ID is unique.

## Ultra-Deep

All FASTQ read IDs will firstly be separated based on their first 2 hex digits of CRC32 and compressed to GZip format.

**WARNING** Extensive memory usage.

## Dependencies

Requires zlib and libdeflate. Also requires `libxxhash-dev` and Intel MKL library.
