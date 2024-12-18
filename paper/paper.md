---
title: 'art_modern: An Accelerated ART Simulator of Diverse Next-Generation Sequencing Reads'
tags:
  - next-generation sequencing
  - illumina sequencer
  - simulation
  - high-performance computing
  - rna-seq
authors:
  - given-names: Zhejian
    surname: Yu
    orcid: 0009-0005-0208-5452
    affiliation: "1"
affiliations:
 - name: The Zhejiang University-University of Edinburgh Institute, Haining (314400), Zhejiang, China
   ror: 04jth1r26
   index: 1
date: 11 December 2024
bibliography: paper.bib
---

# Summary \& Statement of need

High-performance simulation of realistic next-generation sequencing (NGS) reads is a prerequisite for algorithm development and benchmarking. However, most existing simulators are either slow or generate data that does not reflect the real-world error profile of simulators. In a recent benchmark, [Artificial Read Transcription (ART)](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) [@Huang2011] outperforms most of its competitors in aspects of speed and reality of quality scores [@Milhaven2022]. However, ART is long unmaintained, and its optimization techniques are nowadays outdated. Here we introduce [`art_modern`](https://github.com/YU-Zhejian/art_modern), a modern re-implementation of the ART simulator with enhanced performance and functionality. Except all the features of the original ART, our software further supports thread-based parallelism, accelerated random generators, synthesis of BAM files, and support over different sequencing depth for individual contigs.

This software can be used to simulate sequencing data for their own research. Common scenarios include benchmarking of DNA- or RNA-Seq alignment algorithms, test whether the self-built (sc)RNA-Seq pipeline performs good or pressure testing of pipelines on a cluster. This simulator is best suited for GNU/Linux-based [High-End Desktops (HEDTs)](https://www.pcmag.com/encyclopedia/term/hedt) with multiple cores and a fast SSD. However, it may also work on laptops, high-performance clusters (HPCs), or Apple macOS/FreeBSD workstations. We believe with such simulators, the testing and benchmarking of NGS-related bioinformatics algorithms can be largely accelerated.

# Methods

## Simulation mode and input/output

The simulator supports simulation of whole-genome sequencing (WGS), transcriptome sequencing, and template (_i.e._, amplicon) sequencing. Except the first mode that only allows unified coverage across the genome, the simulator supports separated coverage information for each contig on both strands. This allows simulation of differentially expressed genes or isoforms using stranded library construction protocols. The coverage information can be passed through tab-separated values (TSV) formats that contain read ID and coverage, or using the [PBSIM3 Transcripts](https://github.com/yukiteruono/pbsim3/) [@Ono2022] which is a clean TSV format that contains read ID, sequence, and strand coverage information.

Multiple parsers of input files were implemented to accommodate the wide range of formats and simulation strategies acceptable by the simulator. For common FASTA files, the simulator additionally supports [HTSLib](https://github.com/samtools/htslib) [@Bonfield2021] for parsing large genomes with exceptionally large contigs with reasonable memory consumption.

For output, all generated reads together with its alignment information will be stored in an intermediate container class `PairwiseAlignment` (PWA). The instances will then be passed to an output dispatcher, which would trigger registered writers that converts PWAs to their native format (e.g., `std::string` or `bam1_t*`). Each writer internally holds a fixed-sized multi-producer single-consumer [Moody Camel Concurrent Queue](https://github.com/cameron314/concurrentqueue) [@Desrochers2014;@Desrochers2014a] of unique pointers of their native formats. With such, the writers are able to accept PWAs in parallel while writing them in a single thread. The blocking enqueue implemented using a manual try-and-sleep manner. With bulk dequeueing, this implementation is 2 times faster than the queue migrated from Python [@Niksic2020].

Writers of the Sequence Alignment/Map (SAM) and Binary Alignment/Map (BAM) format [@Li2009] are supported through HTSLib, which adds another checking of possibly malformed alignments and parallel compression using arbitrary levels. Except traditional SAM/BAM format, we also support writing of headless SAM/BAM format that does not contain contig information in their heading section with alignment information in the `OA` tag of each alignment record. With such, `art_modern` can perform simulation on a stream of short (c)DNA molecules generated by some high-level simulator. The writer of FASTQ format is implemented using `std::snprintf`, which is considerably faster than `std::stringstream` or Boost.Format. A text-based writer for PWA is also available for debugging purposes.

## Parallelism

To introduce parallelism, a minimal simulation job representation (`SimulationJob` class) was developed. The class contains a job ID, a reference sequence fetcher, and coverage information. Each `SimulationJob` will be passed to an `ArtJobExecutor`, which would generate all reads with the additional parameters it holds. Each thread will hold one `ArtJobExecutor` instance with parallelism realized using the thread pool from Boost.ASIO. A fail-safe mode without any parallelism was also implemented.

Different simulation mode and input parser result in different parallelism strategy. For WGS simulation, each thread gets either one pointer to common in memory reference sequences (`memory` FASTA parser) or have its own pointer to a wrapper of `faidx_t*` (`htslib` FASTA parser), with coverage information being divided and assigned to each thread. For transcriptome and template simulation, however, each thread gets only a part of the reference sequences with the full coverage information (which is small and easier to be duplicated). It is worth noticing that if the reference sequence had already been read into the memory, it will have to be duplicated to be assigned to each thread, which results in 2 times of memory consumption.

## Other Acceleration

Single Instruction, Multiple Data (SIMD) Instruction Set Extensions (ISEs) are widely used acceleration techniques in modern CPUs and C++ Standard Template Libraries (STL). With such technologies, we can achieve a significant performance improvement by grouping similar operations into SIMD-accelerated function calls. For instance, the `for` loops in conversion from reference to read and synthesizing alignment information is accelerated using `std::memcpy` (which is internally accelerated using AVX2 instructions for supporting CPUs in GNU C++ Standard Library) for neighboring matching nucleotides, and successive random number generation is accelerated using bulk generation of random numbers.

A majority of acceleration is done using faster randomization libraries. In this work, we eliminated the usage of `rand` and `srand` function in C Standard Library and replaced them using Mersenne Twister pseudo-random number generator (mt19937) [@Matsumoto1998]. We've also implemented randomization routines using Intel oneAPI Math Kernel Library (oneMKL), which is considerably faster as it supports SIMD-accelerated bulk generation of random numbers. We still kept randomization routines implemented using STL, Boost.Random, and GNU Science Library (GSL) [@Gough2009] for compatibility and portability reasons.

Some functions concerning conversion between raw quality values and their ASCII representation were accelerated using Streaming SIMD Extensions 2 (SSE2) instructions (which should be available in almost all modern Intel CPUs). However, since those functions only run on sequences of limited length and occupied only a small portion of the total CPU time, their acceleration to the entire program is not significant. Further acceleration is possible by using SIMD ISEs with larger vector widths (e.g., AVX2, AVX-512). Other acceleration techniques like copy elimination are also applied.

# Benchmarks

`art_modern` is benchmarked against the original ART 2016.06.05, [DWGSIM](https://github.com/nh13/DWGSIM) 0.1.15, and `wgsim` simulator bundled with [SAMtools](https://github.com/samtools/samtools) [@Li2009;@Danecek2021] 1.21. All tools were compiled using Intel oneAPI C++/DPC++ compiler 2025.0 with `-mtune=native` and `-Ofast` params. `wgsim` and `art_modern` were linked to HTSLib version 1.21 compiled with the same compiler and flag with the original ART linked to GSL 2.8 compiled with the same compiler and flag.

The testing data are the _C. elegans_ reference genome (ce11) from UCSC (hereby referred to as the "genome" data) with a unified sequencing depth of 10, and the long mRNA sequences from UCSC (hereby referred to as the "transcriptome" data) with a unified sequencing depth of 4. DWGSIM were not tested on transcriptome data since it failed with an error `failed to generate a read after 10001 trials`, which may due to the number of ambiguous bases in the reference sequence.

To ensure fair measurement, all simulators were asked to generate reads with PE150 library conformation 3 times on each dataset. The variation generation function of `wgsim` and DWGSIM simulator turned off, and the DWGSIM source code was modified to supress compression of the output FASTQ files. [GNU Time](https://www.gnu.org/software/time/) 1.9 was used to measure the performance of each tool. R [@RCT2024] 4.2.2 with dplyr [@Wickham2023] 1.1.4 and ggplot2 [@Wickham2016] 3.4.4 were used to generate the plots. The benchmark code is provided in the same repository.

From the following figure, we may observe that the CPU time of `art_modern` in both data are significantly lower than that of the original ART and is comparable to DWGSIM, and its wall-clock time is also significantly lower thanks to parallelism. It is worth noting that `art_modern` consumes more memory in the transcriptome data, which may due to its batch-based parallelism strategy while other simulators simulate reads on reference contigs one by one.

![Performance of each tool on two datasets]()

We also compared the generated read qualities using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.12.1, which shows that the distribution of generated read qualities are almost identical to the original ART under the same quality distribution. The quality distribution of `wgsim` and DWGSIM does not reflect real Illumina quality distributions, so is not recommended for real-world applications.

![Quality]()

# Conclusion

In this work, we modified the original ART code base to achieve a significant performance improvement, extensive scalability, and a lot more functions. The resulting software can be served as a drop-in replacement for the original ART. The acceleration strategy and data structures involved can be applied to a wider range of simulators, including long-read simulators like PBSIM3. The software is freely available in [GitHub](https://github.com/YU-Zhejian/art_modern) under the [GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html) license, same as original ART.

# Acknowledgements

The author would like to thank his supervisor, Dr. Wanlu LIU, without whose support this work would not be possible.

# References
