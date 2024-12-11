# ART

Converted from <https://www.niehs.nih.gov/research/resources/software/biostatistics/artusing [pandoc](https://pandoc.org/) with extensive manual modifications.

## Set of Simulation Tools

ART is a set of simulation tools to generate synthetic next-generation
sequencing reads. ART simulates sequencing reads by mimicking real
sequencing process with empirical error models or quality profiles
summarized from large recalibrated sequencing data. ART can also
simulate reads using user own read error model or quality profiles. ART
supports simulation of single-end, paired-end/mate-pair reads of three
major commercial next-generation sequencing platforms: Illumina\'s
Solexa, Roche\'s 454 and Applied Biosystems\' SOLiD. ART can be used to
test or benchmark a variety of method or tools for next-generation
sequencing data analysis, including read alignment, de novo assembly,
SNP and structure variation discovery. ART was used as a primary tool
for the simulation study of the [1000 Genomes
Project](https://www.internationalgenome.org/){target="_blank"}. ART is
implemented in C++ with optimized algorithms and is highly efficient in
read simulation. ART outputs reads in the FASTQ format, and alignments
in the ALN format. ART can also generate alignments in the SAM alignment
or UCSC BED file format. ART can be used together with genome variants
simulators (e.g.
[VarSim](http://bioinform.github.io/varsim/){target="_blank"}) for
evaluating variant calling tools or methods.

## Citation

-   [Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a
    next-generation sequencing read simulator, Bioinformatics (2012) 28
    (4):
    593-594](https://doi.org/10.1093/bioinformatics/btr708)

## Availability

ART is freely available to public. The binary packages of ART are
available for three major operating systems: Linux, Macintosh, and
Windows. ART is also available as Platform-independent C++ source
packages. Each package includes programs, documents and usage examples.

## ChangeLog

### ART-MountRainier-2016-06-05 (the latest version)

This is a minor maintenance release. Thanks to Rene Lujan and other
users for their feedback.

-   fixed a bug leading to the indel rates of 2nd reads same as those
          of 1st reads
-   added a parameter to allow user to set the maximum number of
          insertions and deletions within a read
-   optimized insertion/deletion simulation
-   other minor changes

### ART-GreatSmokyMountains-04-17-2016

This release added the support of the latest Illumina sequencing
systems, fixed a crash bug, and made several other improvements as
detailed below. While ART has always been free to the public, the
release officially puts all ART codes under the [GPL version 3
license](https://www.gnu.org/licenses/gpl-3.0.html){target="_blank"}.
In addition, thanks to Andreas Tille and his Debian Med team, ART is
in the process to be included in the Debian operating system. Thanks
also to many users for their feedback.

-   added built-in quality profiles for the following Illumina
          sequencing systems
          1.  HSXn - HiSeqX PCR free, v2.5, 150bp
    2.  HSXt - HiSeqX TruSeq, v2.5, 150bp
    3.  MinS - MiniSeq TruSeq, 50bp, single-end only
    4.  MSv3 - MiSeq v3, 250bp
    5.  NS50 - NextSeq500 v2, 75bp
-   added the command line option to set the minimum (\--minQ) and
          maximum (\--maxQ) base quality scores for Illumina read simulation
-   fixed a bug related to the crash when simulating Illumina reads
          with a high error rate
-   updated illumina quality profile builder (art_profiler_illumina)
          to make it better and more robust The script now works on both
          MacOS and Linux/Unix systems, and allows users to set the number
          of threads to run the program
-   put all ART codes under [GPL version 3
          license](https://www.gnu.org/licenses/gpl-3.0.html)
-   made several other minor changes and enhancements

### ART-ChocolateCherryCake-03-19-2015

This is a minor maintenance release with two updates thanks to the
feedback from Yancy Lo.

-   art_illumina, art_454, and art_SOLiD
          1.  correct a typo in the HQ tag field of ART SAM alignment files
    2.  add an option to use \'M\' instead of the default \'=/X\' for
                  alignment match/mismatch CIGAR in ART SAM alignment files

### ART-ChocolateCherries-03-09-2015

This release mainly updated ART Illumina simulator with all changes
listed below. Thanks many users for their feedback and contributions
to make the release possible.

-   art_illumina
          1.  fix the issue that longer or shorter reads are produced in the
                        new error-free SAM file when the original reads having indels
    2.  fix the issue that generates empty reads when a reference
                  sequence is too short by skipping the reference
    3.  resolve the problem of doing endless looping when a reference
                  sequence is much shorter than the defined mean fragment length
    4.  fix a bug related to read length in ART_profiler_Illumina
                  (thanks to Thomas Hack)
    5.  add an option \"-c \--rcount\" to specify #number of
                  reads/read pairs to be simulated
    6.  add built-in profiles of HiSeq1000 (100bp), HiSeq 2000(100bp)
                  and HiSeq 2500(125bp and 150bp)
    7.  add an option to specify a sequencing system of which a
                  built-in profile is used for simulation
    8.  report the sequencing system name used for simulation in ART
                  output summary log


### ART-VanillaIceCream-03-11-2014

This is a major update release including several new/enhanced features
and bug fixes. All changes are listed blow. In addition, the accessory
tools ART_profiler_454 and ART_profiler_illumina for ART read profile
generation are included in each distribution package, so it is no need
to download them separately. Thanks all users for their feedback and
contributions that make the release possible.

-   art_illumina
          1.  add the function of amplicon sequencing simulation
    2.  add the support to generate sam file with zero-sequencing
                  errors
    3.  change the maximum allowed quality score to be 93 when scaling
                  quality scores
    4.  fix a bug related to simulation with a fixed random seed for
                  paired-end simulation
    5.  solve the crash issue when running on the newer MacOS X (10.9)
    6.  allow reference sequences to be both DNA and RNA
    7.  fix a bug related to masked \'N\' regions when having multiple
                  reference sequences
    8.  add a cutoff frequency of \"N\"s of genomic regions for
                  masking
    9.  enable turning off the masking of \'N\' regions
    10. other small improvements
-   art_454
          1.  add amplicon sequencing simulation function
    2.  allow reference sequences to be RNA sequences
    3.  enable using a fixed random seed for simulation
    4.  fix the paired-end read direction issue
    5.  enhance log report and other improvements
-   art_SOLiD
          1.  enable simulation with a fixed random seed to generate
                        identical datasets from two simulations
    2.  add support of F3-F5 paired-end read simulation
    3.  add support of amplicon sequencing simulation
    4.  allow reference sequences to be RNA
    5.  enable simulation of reads up to 75bp in length with a testing
                  error profile
    6.  fix the issue of masking \"N\" genomic regions
    7.  make it more user-friendly and some other improvements


### ART-GrapeWine-08-15-2012

The release mainly fixed bugs for ART 454 simulator and added new
features and functions as listed in the following:

1.  support GS FLX Titanium platform
2.  provide new built-in 454 read profiles for both GS FLX and GS FLX
          Titaium platforms
3.  add a new tool 454_readprofile_art that allows users to generate
          their own read profiles from new 454 read data
4.  add an option to change the default flow cycle number
5.  change to not output ALN files by default
6.  change the output of DNA sequences from lower case to upper case
7.  switch automatically to the default indel error profile when user
          own profile does not provide it

### ART-PeachPie-05-16-2012

### ART-CoconutCoffee-01-10-2012

### ART-CranberryJuice-11-23-2011

### ART-ApplePie-04-21-2011
### Installation

**Compilation and installation from a source package**

Compilation of ART from its source codes requires the GNU Scientific
Library (GSL). The GSL can be freely downloaded from GNU at [GSL
Software Website](https://www.gnu.org/software/gsl/){target="_blank"}
. To compile under a Linux/Unix-like operating system, please first
download and unpack a desired source package, then enter the first
directory of the unpacked package, and issue the following commands:

                    ./configure
                    make
                    make install

**Installation from a binary package**

Installation is to simply unpack the binary package to your
installation directory. The executable programs are art_454,
art_illumina, and art_SOLiD for 454, Illumina, and SOLiD platforms,
respectively. Under Linux or MacOS, please use the following command
to unpack a \*.tar.gz binary package:

                    tar xfz art_*.tar.gz

ART binary package for Windows OS is in a ZIP package. You can right
click a ZIP package, and click \"extract\" in the context-menu to
unpack the package.

## Usages

Simple ART usages and examples are given below. Please refer to the
README file in each distribution package for examples and other detail
documentation.

**454 read simulation**

1.  Single-end reads\
    `art_454 [ -s ] [ -p read_profile ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE>`\
    Example:\
    `art_454 seq_reference.fa ./outdir/dat_single_end 20`
2.  Paired-end reads\
    `art_454 [ -s ] [ -p read_profile ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DE>`\
    Example:\
    `art_454 seq_reference.fa ./outdir/dat_paired_end 20 500 20`

**Illumina read simulation**

1.  Single-end reads\
    `art_illumina [options] -i <INPUT_SEQ_FILE> -l <READ_LEN> -f <FOLD_COVERAGE> -o <OUTPUT_FILE_PREFIX>`\
    Example:\
    `art_illumina -sam -i seq_reference.fa -l 50 -f 10 -o ./outdir/dat_single_end`
2.  Paired-end reads\
    `art_illumina [options] -i <INPUT_SEQ_FILE> -l <READ_LEN> -f <FOLD_COVERAGE> -o <OUTPUT_FILE_PREFIX> -m <MEAN_FRAG_LEN> -s <STD_DE>`\
    Example:\
    `art_illumina -p -sam -i seq_reference.fa -l 50 -f 20 -m 200 -s 10 -o d./outdir/dat_paired_end`
3.  Mate-pair reads\
    `art_illumina [options] -i <INPUT_SEQ_FILE> -l <READ_LEN> -f <FOLD_COVERAGE> -o <OUTPUT_FILE_PREFIX> -m <MEAN_FRAG_LEN> -s <STD_DE>`\
    Example:\
    `art_illumina -mp -sam -i seq_reference.fa -l 50 -f 20 -m 2050 -s 50 -o d./outdir/dat_paired_end`

**SOLiD read simulation**

1.  Single-end reads\
    `art_SOLiD [ -s ] [ -p read_profile ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <READ_LEN> <FOLD_COVERAGE>`\
    Example:\
    `art_SOLiD -s seq_reference.fa ./outdir/dat_single_end 32 10`
2.  Paired-end reads\
    `art_SOLiD [ -s ] [ -p read_profile ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <READ_LEN> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DE>`\
    Example:\

    `art_SOLiD seq_reference.fa ./outdir/dat_paired_end 25 10 500 20`

## Contact

Last Reviewed: April 09, 2021
