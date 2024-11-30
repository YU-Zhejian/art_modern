# Frequently Asked Questions (FAQs)

## Supported CPUs?

Although this application should theoretically support both endianness, only little endian is tested. That is, if you're working on an Intel or AMD CPU, this application should work fine.

## How to split produced pair-end/mate-pair sequencing results to 2 FASTQ files?

This can be done through [`seqtk`](https://github.com/lh3/seqtk). For example, to split `tmp/test_small_pe/NC_001416.1.fq`:

```shell
# Read 1
seqtk seq tmp/test_small_pe/NC_001416.1.fq -1 > tmp/test_small_pe/NC_001416.1_1.fq
# Read 2
seqtk seq tmp/test_small_pe/NC_001416.1.fq -2 > tmp/test_small_pe/NC_001416.1_2.fq
```

You may also generate SAM/BAM files and extract PE FASTQ from them using `samtools`.

## Is there an interface for Python, Java, and more?

We may develop a C interface in the future after the APIs and design of the core library are settled.

## How to add adaptors \& primers to the reads?

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adapters to them, and use `templ` mode to introduce sequencing errors.

## I don't care about portability. How to make it wicked fast?

Easy. Set `USE_HTSLIB` to the latest HTSLib available on your system, `CMAKE_BUILD_TYPE` to `Release` or `RelWithDebInfo`, and `USE_RANDOM_GENERATOR` to `ONEMKL` on Intel/AMD machines. Please also make sure that your HTSLib has been linked with [libdeflate](https://github.com/ebiggers/libdeflate).

Also, when executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.
