# Frequently Asked Questions (FAQs)

## Who are the intended user of this simulator?

The intended user of this simulator is anyone who wants to simulate sequencing data for their own research. Typical scenarios includes:

- Benchmarking of DNA- or RNA-Seq alignment algorithms.
- Test whether the RNA-Seq pipeline built by your lab performs well.
- Perform pressure testing of pipelines on a cluster.

This simulator would be best suited for GNU/Linux-based [High-End Desktops (HEDTs)](https://www.pcmag.com/encyclopedia/term/hedt) with multiple cores and a fast SSD. However, it can also work on Laptops, or high-performance clusters (HPCs) with only one node. The current version does not have cross-node parallelism, but we're working on supporting this through MPI.

## What kinds of CPUs are supported by the simulator?

Although this application should theoretically support both endianness, only little endian is tested. That is, if you're working on an Intel or AMD CPU, this application should work fine.

## How to split produced pair-end/mate-pair sequencing results to 2 FASTQ files?

This can be done through [`seqtk`](https://github.com/lh3/seqtk). For example, to split `1.fq`:

```shell
# Read 1
seqtk seq 1.fq -1 > 1_1.fq
# Read 2
seqtk seq 1.fq -2 > 1_2.fq
```

You may also generate SAM/BAM files and extract PE FASTQ from them using `samtools`.

## Is there an interface for Python, Java, and more?

We may develop a C interface in the future after the APIs and design of the core library are settled.

## How to add adaptors \& primers to the reads?

Currently, there's no support for such features in the simulator. However, you may manually chop your reference genome, add adapters to them, and use `templ` mode to introduce sequencing errors.

## My Question is Still Unanswered

Please feel free to open an issue on GitHub.
