---
title: 'art_modern: A Modern Simulator of Diverse Next-Generation Sequencing Reads'
tags:
  - next-generation sequencing
  - illumina sequencer
  - simulation
  - TODO1 # TODO
  - TODO2 # TODO
authors:
  - given-names: Zhejian
    surname: Yu
    orcid: 0000-0000-0000-0000
    affiliation: "1"
affiliations:
 - name: The Zhejiang University-University of Edinburgh Institute, Haining (314400), Zhejiang, China
   index: 1
date: 11 December 2024
bibliography: paper.bib
---

# Summary \& Statement of need

High-performance simulation of realistic next-generation sequencing (NGS) data is a must for various algorithm development and benchmarking tasks. However, most existing simulators are either slow or generate data that does not reflect the real-world error profile of simulators. Here we introduces [`art_modern`](https://github.com/YU-Zhejian/art_modern), a modern re-implementation of the popular [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) [@Huang2011] simulator with enhanced performance and functionality. Except all the features of the original ART, our software further supports thread-based parallelism, accelerated random generators, synthesis of BAM files, and support over different sequencing depth for individual contigs.

This software can be used to simulate sequencing data for their own research. Common scenarios include benchmarking of DNA- or RNA-Seq alignment algorithms, test whether the self-built RNA-Seq pipeline performs well or pressure testing of pipelines on a cluster. This simulator is best suited for GNU/Linux-based [High-End Desktops (HEDTs)](https://www.pcmag.com/encyclopedia/term/hedt) with multiple cores and a fast SSD. However, may can also work on Laptops, high-performance clusters (HPCs) with only one node, or Apple macOS and FreeBSD workstations. We believe with such simulators, the testing and benchmarking of NGS-related bioinformatics algorithms can be largely accelerated.

# Methods

## Upper-Level Design

## I/O

The SAM and BAM output format is supported through [HTSLib](https://github.com/samtools/htslib) [@Bonfield2021].

The asynchronous output writer is implemented using [Moody Camel Concurrent Queue](https://github.com/cameron314/concurrentqueue) [@Desrochers2014].

# Benchmarks



# Availability

The software is available in [GitHub](https://github.com/YU-Zhejian/art_modern) under the [GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html) license, same as original ART.

# Acknowledgements

The author would like to thank his supervisor, Dr. Wanlu LIU, without whose support this work would not be possible.

# References
