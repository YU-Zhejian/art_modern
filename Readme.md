# `art_modern`: Modernized ART Simulator of Diverse Next-Generation Sequencing Reads

[![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/YU-Zhejian/art_modern/total.svg)](https://github.com/YU-Zhejian/art_modern/releases/)
[![License](https://img.shields.io/badge/licence-GPL_3.0-blue.svg)](https://www.gnu.org/licenses/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/art_modern/README.html)
[![Anaconda.org](https://anaconda.org/bioconda/art_modern/badges/version.svg)](https://anaconda.org/bioconda/art_modern)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/art_modern/badges/downloads.svg)](https://anaconda.org/bioconda/art_modern)

[![Catalogued on GSR](https://surveillance.cancer.gov/gsr/static/img/gsr_tile.jpg)](https://surveillance.cancer.gov/genetic-simulation-resources/packages/art_modern/)

## Introduction

Here we introduce `art_modern`, a modern re-implementation of the popular [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) simulator with enhanced performance and functionality. It can be used for anyone who wants to simulate sequencing data for their own research, like benchmarking of DNA- or RNA-Seq alignment algorithms, test whether the RNA-Seq pipeline built by your lab performs well or perform pressure testing of pipelines on a cluster.

## Quick Start

### Installation

**NOTE** Although we're providing compiled binaries for Linux, we still recommend compiling the source code yourself for maximum performance.

#### Using Fully Static Build

In each release, there will be a file named `build_rel_with_dbg_alpine-x86_64.zip`/`build_rel_with_dbg_alpine-x86_64.tar.gz` in the [Releases](https://github.com/YU-Zhejian/art_modern/releases) section. The file contains fully static linked libraries and executable binaries built under x86\_64 Alpine Linux, that should work on most x86\_64 Linux distributions. Unzip it and you're good to go.

**WARNING** Static builds may lead to compromized security.

#### Using `dpkg`

If you use specific versions of Debian or Ubuntu, you can install `art_modern` through pre-built DEB packages available with each release. Root privileges are required.

#### Using Docker/Podman

See: <https://quay.io/repository/biocontainers/art_modern>.

#### Using Conda

[Conda](https://docs.conda.io/) (or [Mamba](https://mamba.readthedocs.io/en/latest/)/[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)) is a popular open-source package and environment management system that simplifies the installation and management of software packages and their dependencies. Before processing, make sure you've installed Conda >=25.7.0 by `conda --version`. Then:

```shell
conda create -y -n art_modern_bioconda -c bioconda -c conda-forge art_modern
```

to create an environment named `art_modern_bioconda` with the package installed.

#### Compiling the Source Code

Use [Git](https://git-scm.com/) to clone this repository:

```shell
git clone https://github.com/YU-Zhejian/art_modern.git
cd art_modern
```

Or, if you have no Git:

```shell
curl -L https://github.com/YU-Zhejian/art_modern/archive/refs/heads/master.zip -o art_modern-master.zip
unzip art_modern-master.zip
cd art_modern-master
```

Ensure you have a C++ compiler that supports [C++17](https://en.cppreference.com/w/cpp/17) installed on your computer (Preferably [GCC](https://gcc.gnu.org/) >= 9.5.0). Also check whether your [CMake](https://cmake.org/), [GNU Make](https://www.gnu.org/software/make/), [Boost C++ Library](https://www.boost.org/), [GNU BinUtils](https://www.gnu.org/software/binutils/), [GNU Bash](https://www.gnu.org/software/bash/), [GNU CoreUtils](https://www.gnu.org/software/coreutils/), [Python 3](https://www.python.org/), and minimal HTSLib-dependencies (namely, [zlib](https://www.zlib.net/) and [pthread](https://www.man7.org/linux/man-pages/man7/pthreads.7.html)) are working.

Build the project using:

```shell
mkdir -p opt/build_release
env -C opt/build_release cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo "$(pwd)"
cmake --build opt/build_release -j40
```

The project binary will be available at `opt/build_release/art_modern`. Now we can test whether the program runs:

```shell
opt/build_release/art_modern --help
opt/build_release/art_modern --version
```

### Simulating WGS Data using _E. Coli_ Genome

Download _E. coli_ reference genome from NCBI. Here we'll use K12 strand MG1655 sub-strand as an example.

```shell
mkdir -p opt/data/
wget \
    -4 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz \
    -O opt/data/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip -k opt/data/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Now we can simulate WGS data using E. coli reference genome. Let's start with single-end sequencing using HiSeq 2500 with 125 bp read length and 10X coverage.

```shell
opt/build_release/art_modern \
    --mode wgs \
    --lc se \
    --i-file opt/data/GCF_000005845.2_ASM584v2_genomic.fna \
    --o-fastq opt/build_release/e_coli_wgs_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov 10
```

The generated FASTQ file will be at `opt/build_release/e_coli_wgs_se.fastq`.

We may also simulate paired-end data with the following configuration:

```shell
opt/build_release/art_modern \
    --mode wgs \
    --lc pe \
    --i-file opt/data/GCF_000005845.2_ASM584v2_genomic.fna \
    --o-fastq opt/build_release/e_coli_wgs_pe.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov 10 \
    --pe_frag_dist_mean 300 \
    --pe_frag_dist_std_dev 50
```

Please note that we have additionally specified a quality file for read 2 with the mean and standard deviation of fragment lengths.

### Simulating RNA-Seq Data using _C. Elegans_ Transcriptome

Simulating transcriptome is a little more complicated since each cDNA molecule has different counts. Strand-specific library technologies also generate RNA-Seq data on one strand only. You're recommended to use [YASIM](https://pypi.org/project/YASIM) or other high-level simulators to generate expression for each cDNA molecule. You may also easily convert outputs from [featureCounts](https://subread.sourceforge.net/featureCounts.html), [htseq-count](https://htseq.readthedocs.io/en/latest/), [Salmon](https://salmon.readthedocs.io/en/latest/), [Kallisto](https://pachterlab.github.io/kallisto/) or [STAR](https://github.com/alexdobin/STAR) to the format supported by `art_modern`. The unified coverage model (i.e., like WGS) is also supported.

Please note that cDNAs with insufficient length will be ignored. We also do not support circular RNA simulation.

#### Unified Coverage

Following example retries the first 1000 transcripts from reference _C. Elegans_ transcriptome from [UCSC Genome Browser](https://genome.ucsc.edu/) and performs a simulation using 10X unified coverage. You need to install [seqtk](https://github.com/lh3/seqtk) to run this example:

```shell
curl https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/mrna.fa.gz | \
    gzip -cdf | \
    seqtk sample /dev/stdin 1000 > opt/build_release/ce11_mrna_1000.fa

opt/build_release/art_modern \
    --mode trans \
    --lc se \
    --i-file opt/build_release/ce11_mrna_1000.fa \
    --o-fastq opt/build_release/c_elegans_trans_unified_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov 10
```

#### Unstranded Coverage

To simulate data with unstranded coverage information (i.e., same coverage on both strands), you need to provide an additional TSV file with one column of transcript ID and another column of coverage (in floating points). Please note that lines started by `#` will be ignored. An example of the coverage file:

```text
NM_069135	6.695683025425357
NR_056112	5.19437291612395
NR_051843	3.4504965075273137
NR_066512	4.73632003156384
```

The following example generates a coverage file using [GNU AWK](https://www.gnu.org/software/gawk/) with random coverage ranging between 0 and 10 assigned to each cDNA molecule.

```shell
samtools faidx opt/build_release/ce11_mrna_1000.fa
awk 'BEGIN{print "#ID\tCOV";}{printf "%s\t%f\n", $1, (rand()*10);}' \
    < opt/build_release/ce11_mrna_1000.fa.fai \
    > opt/build_release/ce11_mrna_1000.fa.unstranded_cov.tsv

opt/build_release/art_modern \
    --mode trans \
    --lc se \
    --i-file opt/build_release/ce11_mrna_1000.fa \
    --o-fastq opt/build_release/c_elegans_trans_unstranded_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov opt/build_release/ce11_mrna_1000.fa.unstranded_cov.tsv
```

#### Stranded Coverage

To simulate data with stranded coverage information (i.e., coverage on one strand is different from the other), you need to provide an additional TSV file with one column of transcript ID and two other columns of coverage in positive and negative strand (in floating points). An example of the coverage file:

```text
NM_069135	2.3137902802960717	4.381892745129285
NR_056112	3.47140212944225	1.7229707866816995
NR_051843	1.3540475385633155	2.0964489689639985
NR_066512	3.0468993830563917	1.689420648507448
```

Code example:

```shell
awk 'BEGIN{print "#ID\tCOV_POS\tCOV_NEG";}{printf "%s\t%f\t%f\n", $1, (rand()*5), (rand()*5);}' \
    < opt/build_release/ce11_mrna_1000.fa.fai \
    > opt/build_release/ce11_mrna_1000.fa.stranded_cov.tsv

opt/build_release/art_modern \
    --mode trans \
    --lc se \
    --i-file opt/build_release/ce11_mrna_1000.fa \
    --o-fastq opt/build_release/c_elegans_trans_stranded_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov opt/build_release/ce11_mrna_1000.fa.stranded_cov.tsv
```

#### The PBSIM3 Transcripts Input Format

The PBSIM3 Transcripts input format is a 4-column tab-delimited text file with transcript ID, sequence, and coverage on both strands. This file includes both sequence and coverage, so no additional coverage parameter is required. Similarly, sequences with insufficient length and lines started with `#` will be ignored. An example of the transcript input file (Sequences represented as `aaaa`):

```text
NR_056112	3.47140212944225	1.7229707866816995	aaaa
NR_051843	1.3540475385633155	2.0964489689639985	aaaa
NR_066512	3.0468993830563917	1.689420648507448	aaaa
NM_061905	0.9618664937744315	1.3989801728399471	aaaa
NR_054174	3.591258844822635	4.92434801892288	aaaa
```

The following example converts the FASTA file to the PBSIM3 Transcripts input format with the help of [seqkit](https://bioinf.shenwei.me/seqkit) with random coverage generated using GNU AWK.

```shell
seqkit fx2tab opt/build_release/ce11_mrna_1000.fa | \
    awk 'BEGIN{print "#ID\tCOV_POS\tCOV_NEG\tSEQ";}{printf "%s\t%f\t%f\t%s\n", $1, (rand()*5), (rand()*5), $3;}' \
    > opt/build_release/ce11_mrna_1000.fa.pbsim3_trans.tsv

opt/build_release/art_modern \
    --mode trans \
    --lc se \
    --i-file opt/build_release/ce11_mrna_1000.fa.pbsim3_trans.tsv \
    --o-fastq opt/build_release/c_elegans_trans_pbsim3_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-type pbsim3_transcripts
```

### Template-Based Simulation

Template-based simulation is often used to introduce Illumina-specific errors to cDNA molecules generated from some upstream simulator like [CAMPAREE](https://camparee.readthedocs.io/en/latest/). In this mode, single-end reads will be started from the first base of the template while paired-end/mate-pair reads will span the entire template. The template-based simulation mode also supports the PBSIM3 Transcripts format. For example:

```shell
opt/build_release/art_modern \
   --mode template \
   --lc pe \
   --i-file opt/build_release/ce11_mrna_1000.fa.pbsim3_trans.tsv \
   --o-fastq opt/build_release/c_elegans_template_pbsim3_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
   --read_len 125 \
   --parallel 4 \
   --i-type pbsim3_transcripts
```

Please note that the mean and standard deviation of fragment length is not specified since in template-based simulation, a template is considered a fragment.

### Using UNIX Pipelines

With UNIX pipelines, we can redirect the input and output of `art_modern` from/to another file/processes. Following example reads FASTA reference from `/dev/stdin` (Standard Input), and writes compressed FASTQ, PWA, and sorted BAM files.

This example requires [gzip](https://www.gnu.org/software/gzip/), [pigz](https://zlib.net/pigz/), [SAMtools](https://github.com/samtools/samtools), and [xz-utils](https://tukaani.org/xz/).

```shell
zcat opt/data/GCF_000005845.2_ASM584v2_genomic.fna.gz | \
    opt/build_release/art_modern \
    --mode wgs \
    --lc se \
    --i-file /dev/stdin \
    --i-type fasta \
    --i-parser memory \
    --o-fastq >(pigz -p8 -9 -v -cf - > opt/build_release/e_coli_wgs_se.fastq.gz) \
    --o-pwa >(xz -9 -T5 -vv -cf - > opt/build_release/e_coli_wgs_se.pwa.xz) \
    --o-sam >(samtools sort -@9 --write-index -o opt/build_release/e_coli_wgs_se.sorted.bam) \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov 5
```

Please wait for a while for the compression to finish.

## Using MPI

If you're working on a cluster with [MPI](https://www.mpi-forum.org/) installed, you can use the MPI version of `art_modern` to speed up the simulation by executing `art_modern` on different nodes in parallel.

**NOTE** Please consult your cluster administrator/cluster documentation before using MPI version, since some clusters may have special configurations and/or quota limitations.

See also:

- [Documentation for MPI Integration in Slurm](https://slurm.schedmd.com/mpi_guide.html) if your cluster uses [Slurm Workload Manager](https://slurm.schedmd.com/).

Documentations from MPI vendors:

- [Open MPI Documentation](https://www.open-mpi.org/doc/), esecially [Launching MPI applications](https://docs.open-mpi.org/en/v5.0.x/launching-apps/index.html) section.
- [MPICH Documentation](https://www.mpich.org/documentation/).
- [Intel MPI Documentation](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html?wapkw=MPI).

You may also find the following information helpful:

- [OpenHPC Community](https://openhpc.community/) if you'd like to build your own HPC cluster.

### Building MPI-Enabled `art_modern`

Load the MPI environment on your cluster. Then execute:

```shell
mkdir -p opt/build_release-mpi
env -C opt/build_release-mpi cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DWITH_MPI=ON "$(pwd)"
cmake --build opt/build_release-mpi -j40
```

Test whether MPI runtime library is correctly linked to the executable. For example:

```shell
ldd opt/build_release-mpi/art_modern-mpi | grep mpi
# [...]
#         libmpi_cxx.so.40 => /lib/x86_64-linux-gnu/libmpi_cxx.so.40 (0x00007ee9f64b0000)
#         libmpi.so.40 => /lib/x86_64-linux-gnu/libmpi.so.40 (0x00007ee9f637c000
```

indicates that the MPI runtime library are correctly linked. Note that the executable name is now suffixed with `-mpi`.

Test whether the MPI version works:

```shell
mpiexec -n 1 opt/build_release-mpi/art_modern-mpi --version
```

You should be able to see the output MPI standard version and vendor-specific information. For example:

```text
MPI:
        Standard Version: 3.1
        Library Version: Open MPI v4.1.6, [...]
```

### Output Files

Let us run a WGS simulation with 4 MPI processes:

```shell
mpiexec -n 4 opt/build_release-mpi/art_modern-mpi \
    --mode wgs \
    --lc se \
    --i-file opt/data/GCF_000005845.2_ASM584v2_genomic.fna \
    --o-fastq opt/build_release-mpi/e_coli_wgs_se.fastq \
    --builtin_qual_file HiSeq2500_125bp \
    --read_len 125 \
    --parallel 4 \
    --i-fcov 20
```

4 files, namely `e_coli_wgs_se.0.fastq`, `e_coli_wgs_se.1.fastq`, `e_coli_wgs_se.2.fastq`, and `e_coli_wgs_se.3.fastq` will be generated in `opt/build_release-mpi/` directory. Each file corresponds to the output from each MPI process (rank 0 to rank 3). This behavior is the same for other output files (e.g., PWA, SAM/BAM).

**NOTE** In the avove example, the actual number of threads used in simulation will be 16 (4 MPI processes * 4 threads per process) since we've specified `--parallel 4`.

### Other Differences Between MPI and Non-MPI Versions

- The simulator will not be able to support UNIX devices and/or redirctions for input/output files. So, even if your input is enormous, you have to write it to a physical file first. However, as most HPC clusters use distributed file systems, this should not be a big problem.
- Logging issues: Only log messages from rank 0 process will be printed to standard error. Other ranks' log messages will be written to disk (If environment variable `ART_NO_LOG_DIR` is not set) or discarded (If environment variable `ART_NO_LOG_DIR` is set).

## What's Next?

The project provides diverse documentations to satisfy your needs.

- If you want to build the software with different options, see [Install](docs/Install.md).
- For a detailed guide on parameters and their combinations, see [Usage](docs/Usage.md) and [FAQ](docs/FAQ.md).
- For developers, please refer to:
  - [Contributing](docs/Contributing.md) for software engineering tasks and contribution guidelines. See also [Code of Conduct](docs/CODE_OF_CONDUCT.md).
  - [Design](docs/Design.md) for the latest design of the software.
  - [Copying](docs/Copying.md) for third-party libraries and codes used in this project.
  - [News](docs/News.md) for changes over the project.
- For a comparison of this project with other simulators, see [Benchmark](explore/benchmark_other_simulators).
- The original ART documentation is copied to [here](docs/original_art_docs).

## Acknowledgements

This simulator is based on the works of [Weichun Huang](mailto:whduke@gmail.com) _et al._, under [GNU GPL v3](https://www.gnu.org/licenses/) license. The software is originally distributed [here](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) with the following reference:

- W. Huang, L. Li, J. R. Myers, and G. T. Marth, _ART: a next-generation sequencing read simulator_, Bioinformatics (Oxford, England), vol. 28, no. 4, pp. 593--594, Feb. 2012, doi: [10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708).

The bundled HTSLib library used MIT License with the following reference:

- J. K. Bonfield et al., _HTSlib: C library for reading/writing high-throughput sequencing data_, GigaScience, vol. 10, no. 2, p. giab007, Jan. 2021, doi: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007).

Other libraries used in this project are distributed under their own licenses. See [Copying](docs/Copying.md) for details.
