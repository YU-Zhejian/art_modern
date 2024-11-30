# `art_modern`: Simulator of Diverse Next-Generation Sequencing Reads

## Quick Start

Build the project using:

```shell
mkdir -p build_release
env -C build_release cmake -DCMAKE_BUILD_TYPE=Release -DCEU_CM_SHOULD_ENABLE_TEST=FALSE ..
env -C build_release make -j40
```

The project binary will be available at `build_release/art_modern`. Now we can test whether the program runs:

```shell
build_release/art_modern --help
build_release/art_modern --version # For version information
```

### Simulating WGS Data using _E. Coli_ Genome

Download _E. Coli_ reference genome from NCBI. Here we'll use K12 strand MG1655 substrand as an example.

```shell
mkdir -p tutorial_data
wget -4 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Now we can simulate WGS data using E. Coli reference genome. Let's satrt with single-end sequencing using HiSeq 2500 with 125bp read length and 10X coverage.

```shell
build_release/art_modern \
   --mode wgs \
   --lc se \
   --i-file tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna \
   --o-fastq tutorial_data/e_coli_wgs_se.fastq \
   --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
   --read_len 125 \
   --parallel 4 \
   --i-fcov 10
```

The generated FASTQ file will be at `tutorial_data/e_coli_wgs_se.fastq`.

We may also simulate paired-end data with following configuration:

```shell
build_release/art_modern \
   --mode wgs \
   --lc pe \
   --i-file tutorial_data/GCF_000005845.2_ASM584v2_genomic.fna \
   --o-fastq tutorial_data/e_coli_wgs_pe.fastq \
   --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt \
   --qual_file_2 data/Illumina_profiles/HiSeq2500L125R2.txt \
   --read_len 125 \
   --parallel 4 \
   --i-fcov 10 \
   --pe_frag_dist_mean 300 \
   --pe_frag_dist_std_dev 50
```

Please note that we have additionally specified quality file for read 2 with the mean and standard diviation of fragment lengths.

### Simulating RNA-Seq Data using _C. Elegans_ Transcriptome

### Simulating Targeted Amplification Data

### Simulating WGS Data on Large Genomes

## What's Next?

The `art_modern` project provides diverse documentations to satisfy your needs.

- If you want to build the software with different options, see [Install](docs/Install.md).
- For detailed guide on parameters and their combinations, see [Usage](docs/Usage.md).
- For developers, please refer to:
  - [Contributing](docs/Contributing.md) for design principle and contribution guidelines.
  - [Copying](docs/Copying.md) for third-party libraries and codes used in this project.
  - [News](docs/News.md) for changes over the project.

## Acknowledgements

This simulator is based on the works of [Weichun Huang](mailto:whduke@gmail.com) _et al._, under [GNU GPL v3](https://www.gnu.org/licenses/) license. The software is originally distributed [here](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) with the following reference:

- W. Huang, L. Li, J. R. Myers, and G. T. Marth, _ART: a next-generation sequencing read simulator_, Bioinformatics (Oxford, England), vol. 28, no. 4, pp. 593–594, Feb. 2012, doi: [10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708).

The bundled HTSLib library used MIT License with the following reference:

- J. K. Bonfield et al., _HTSlib: C library for reading/writing high-throughput sequencing data_, GigaScience, vol. 10, no. 2, p. giab007, Jan. 2021, doi: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007).
