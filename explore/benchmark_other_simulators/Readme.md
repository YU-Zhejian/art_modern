# Benchmarks

Benchmark of `art_modern` with WGSim (Bundled with SAMtools) and the original ART simulator. All benchmarks are performed on a Linux machine with Intel DPC++/C++ compiler.

## Instructions

Use `build.sh` to build other simulators.

Use `build-update.sh` to build and update `art_modern`.

Use `run.sh` to run benchmark.

## Datasets

1. Reference genome of _C. Elegans_.
2. Long mRNAs of human that exceeds 500 nucltotides.

The simulators were asked to generate a 10X sequencing data of FASTQ format on both datasets with paired-end reads lengthen 150bp.

## Results

### PC 1

- **MODEL** HP ENVY x360 Convertible 15 dr1xxx
- **CPU** 13th Gen Intel(R) Core(TM) i7-13700H
- **MEM** TODO
- **SSD** TODO
- **OS** Linux Mint 22
- **KERNEL** Linux 6.8.0-49-generic (x86_64) #49-Ubuntu SMP PREEMPT_DYNAMIC Mon Nov 4 02:06:24 UTC 2024

### PC 2

- **MODEL** Self-assembled HEDT.
- **CPU** TODO
- **MEM** TODO
- **SSD** TODO
- **OS** TODO
- **KERNEL** TODO

## Used Third-Party Codes

### SAMtools 1.21 Release Tarball

Available from <https://sourceforge.net/projects/samtools/files/samtools/1.21/samtools-1.21.tar.bz2>

Affected files:

- `src/wgsim.c`

With MIT/Expat License.

### Original ART Source Code v2016.06.05

Available from <https://www.niehs.nih.gov/sites/default/files/2024-02/artsrcmountrainier2016.06.05linux.tgz>

Affected files:

- `src/art_original/*`

With GNU GPL v3 License.
