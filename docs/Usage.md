# Usage

This is the detailed usage of `art_modern`. Every parameter and their combinations are introduced in detail. In this documentation, commands will be represented as `ls -lFh` with in-line or block omission represented as `[...]`.

## Simulation Modes (`--mode`)

- **`wgs` (DEFAULT) for whole-genome sequencing.** For scenarios with a limited number of sized contigs.
- `trans` for transcriptome simulation. For scenarios with many short contigs.
- `template` for template simulation. Often being used with a high-level simulator that produces read fragments.

## Library Construction Methods (`--lc`)

See [here](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) and [here](https://www.illumina.com/science/technology/next-generation-sequencing/mate-pair-sequencing.html) for details on library construction methods.

- **`se` (DEFAULT) for single-end reads.** That includes:
  - **For `wgs` and `trans` mode:** A read may start from any position in the reference contig that is long enough to contain the read.
  - **For `template` mode:** A read will start from 0 position in the reference contig.
- `pe` for paired-end reads. That includes:
  - **For `wgs` and `trans` mode:** A fragment, whose length is drawn from a Gaussian distribution specified by `--pe_frag_dist_mean` and `--pe_frag_dist_std_dev`, is created from any position in the reference contig that is long enough to contain the fragment.
  - **For `template` mode:** A fragment that spans the entire length of the reference contig is created.
  - Paired-end reads facing inward of the fragment are then created from both ends of the fragment.
- `mp` for mate-paired reads. That includes:
  - Extraction of the fragment as-is in `pe` mode.
  - Mate-paired reads facing outward of the fragment are created from both ends of the fragment.

Parameters `--pe_frag_dist_mean` and `--pe_frag_dist_std_dev` are needed to specify the length distribution of fragments. In original ART, the fragment lengths obey Gaussian (normal) distribution.

## Parallelism (`--parallel`)

Number of threads to use. Use `-1` to disable parallelism. Use `0` (default) to use all available cores.

## Input (`--i-*`)

Currently, we support input in FASTA and PBSIM3 Transcripts format. They are controlled by the following major parameters:

- `--i-file`: The input reference file path. It must exist in the filesystem. Device paths like `/dev/stdin` is supported.
- `--i-type`: The file type of input reference sequences. Currently, FASTA and PBSIM3 Transcripts format are supported.
  - **`auto` (DEFAULT) for extension-based decision.**
    - If the file ends with `.fna`, `.fsa`, `.fa`, `.fasta`, resolve to `fasta`.
    - Otherwise, an error will be raised.
  - `fasta` for FASTA files.
  - `pbsim3_transcripts` for PBSIM3 Transcripts format.

### Input Parser (`--i-parser`)

- **`auto` (DEFAULT) for size-based determination.**
  - If the file size larger than 1 GiB or can not be told (which is quite common if the input was redirected from stdin or other devices), resolve to `htslib` (`wgs` mode) or `stream` (`trans` or `template` mode).
  - Otherwise, use `memory`.
- `memory`: The entire file will be read into the memory.
  - Fast for small reference files.
- `htslib`: Store FASTA Index in memory and fetch sequences from disk using [HTSLib](https://github.com/samtools/htslib).
  - A FASTA Index (Usually ends with `.fai` and can be built using `samtools faidx`) will be needed.
  - Memory-efficient for large genome assembly FASTA files with limited number of long contigs.
  - Inefficient for transcriptome/template FASTAs with large number of (relatively) short contigs.
  - Each thread will hold its own FASTA Index in memory.
  - **DO NOT SUPPORT PBSIM3 TRANSCRIPT FORMAT.**
- `stream`: Streamline the input reference as batches and process them one by one.
  - Efficient for transcriptome/template FASTAs with large number of (relatively) short contigs.
  - **CONTIG NUMBER AND LENGTH INFORMATION NOT AVAILABLE**, so use headless SAM/BAM output if a SAM/BAM output is needed.
  - Batch size controlled by `--i-batch_size`.

### More Instructions on FASTA Format

FASTA format can be parsed by all parsers. However, please keep in mind that:

**NOTE** For `htslib` parser, identical line lengths (except the last line) inside a contig is assumed. That means the following FASTA file is legal for `htslib` parser:

```text
>chr1
AAAAAAAAAAAA
AAAA
>chr2
AAAAAAAAAAAA
AAAAAAAAAAAA
AAAAAAAAAAAA
AA
>chr3
AA
```

But the following is not:

```text
>chr2
AAAA
AAAAAAAAAA
AAAAAAAA
AA
```

**NOTE** For read names, only characters before the first whitespace characters (space ` `, tabs `\t`, etc.) are read. That is, the FASTA file:

```text
>chr1 some attrs
AAAAAATTTTTT
>chr2 more attrs
AAAAAATTTTTT
```

Will be parsed into identical data structure with:

```text
>chr1
AAAAAATTTTTT
>chr2
AAAAAATTTTTT
```

### Coverage (`--i-fcov`)

**NOTE** This parameter is ignored if the file type is set to `pbsim3_transcripts`.

This option allows you to specify coverage. `art_modern` supports the following coverage mode:

- Unified coverage in `double` data type (e.g., 10.0). Under this scenario, the coverage is identical for all contigs.
- Per-contig coverage without strand information.
- Per-contig coverage with strand information.

**NOTE** We prefer unified coverage. That is, if you set this parameter to 10.0, we're **NOT** going to check whether there's a file named `10.0`. So please make sure that the file name is not a number if you want to specify per-contig coverage.

### Conclusive Remarks

Compatibility matrix of file type, simulation mode, and parser:

| Parser \ Simulation Mode | `wgs`     | `trans`                        | `template`                     |
|--------------------------|-----------|--------------------------------|--------------------------------|
| `memory`                 | `fasta`   | `fasta` / `pbsim3_transcripts` | `fasta` / `pbsim3_transcripts` |
| `htslib`                 | `fasta`   | **ERROR**                      | **ERROR**                      |
| `stream`                 | **ERROR** | `fasta` / `pbsim3_transcripts` | `fasta` / `pbsim3_transcripts` |

Compatibility matrix of coverage mode, simulation mode, and file type:

| Simulation Mode \ File Type | `fasta`                        | `pbsim3_transcripts` |
|-----------------------------|--------------------------------|----------------------|
| `wgs`                       | Unified                        | **ERROR**            |
| `trans`                     | Unified / Unstraded / Stranded | **IGNORED**          |
| `template`                  | Unified / Unstraded / Stranded | **IGNORED**          |

## Output Formats (`--o-*`)

Here introduces diverse output formats supported by `art_modern`. You may specify none of them to perform simulation without any output for benchmarking purposes.

If the parent directory does not exist, it will be created using a `mkdir -p`-like manner.

### Pairwise Alignment Format (`--o-pwa`)

PWA is a serialization of the internally used data structure of the simulator. Currently, PWA consists of some metadata lines (Started with `#`) and a list of alignments, each spanning 4 lines.

The alignment lines are as follows:

1. `>`, read name, a blank, and the leftmost position of the alignment.
2. The aligned read sequence, gapped with `-`.
3. The aligned reference sequence, gapped with `-`.
4. The gapless quality string.

Example:

```text
#PWA
#ARGS: opt/build_debug/art_modern --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt [...]
>NM_069135:art_modern:1:nompi:0	NM_069135:0:+
AGC-CAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGT-TCGAA
AGCACAAA-GGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCTTGAGAGGTAGTTAGAGAAAACGAGAAA-ACTACGAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTTCGAA
BCCCCGGGGGGGFGGGGGGGFGGGGGGG1GGGGGGGFG1GGFGGG:GGG/GDGGGGGGG:GGGEGGGGGCGGGGGGGGEGGD<DGGGGGGDF>GGGG0GG:FGGGGGGGGGGCG.GEGGGGGGGG
[...]
```

### FASTQ Format (`--o-fastq`)

The good old FASTQ format. The qualities are Phread encoded in ASCII with an offset of 33.

```text
@NM_069135:art_modern:1:nompi:0
AGCCAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTCGAA
+
BCCCCGGGGGGGFGGGGGGGFGGGGGGG1GGGGGGGFG1GGFGGG:GGG/GDGGGGGGG:GGGEGGGGGCGGGGGGGGEGGD<DGGGGGGDF>GGGG0GG:FGGGGGGGGGGCG.GEGGGGGGGG
```

FASTQ files can be easily converted to FASTA using [`seqtk`](https://github.com/lh3/seqtk):

```shell
cat in.fq | seqtk seq -A > out.fa
```

See also:

- [Specifications of Common File Formats Used by the ENCODE Consortium at UCSC](https://genome.ucsc.edu/ENCODE/fileFormats.html#FASTQ).
- [Common File Formats Used by the ENCODE Consortium](https://www.encodeproject.org/help/file-formats/#fastq)

### SAM/BAM Format (`--o-sam`)

Sequence Alignment/Map (SAM) and Binary Alignment/Map (BAM) format supports storing of ground-truth alignment information and other miscellaneous parameters. They can be parsed using [samtools](https://github.com/samtools/samtools), [`pysam`](https://pysam.readthedocs.io/) and other libraries, and can be used as ground-truth when benchmarking sequence aligners.

This writes canonical SAM/BAM format using [HTSLib](https://github.com/samtools/htslib). This output writer supports computing `NM` and `MD` tag. The mapping qualities for all reads are set to 255.

Example:

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt [...]
@SQ	SN:NM_069135	LN:1718
[...]
NM_069135:art_modern:1:nompi:0	0	NM_069135	1	255	3=1D4=1I41=1X27=1I5=1X36=1D5=	=	1	0	AGCCAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTCGAA	cddddhhhhhhhghhhhhhhghhhhhhhRhhhhhhhghRhhghhh[hhhPhehhhhhhh[hhhfhhhhhdhhhhhhhhfhhe]ehhhhhheg_hhhhQhh[ghhhhhhhhhhdhOhfhhhhhhhh	MD:Z:3^A45T32G36^T5	NM:i:6
[...]
```

Other SAM/BAM formatting parameters includes:

- `--o-sam-use_m`: Computing CIGAR string using `M` (`BAM_CMATCH`) instead of `=`/`X` (`BAM_CEQUAL`/`BAM_CDIFF`). Rarely used in next-generation sequencing but common for long-read sequencing alignments.
- `--o-sam-write_bam`: Write BAM instead of SAM.
- `--o-sam-num_threads`: Number of threads used to compress BAM output.
- `--o-sam-compress_level`: [`zlib`](https://www.zlib.net/) compression level. Supports `[u0-9]` with `u` for uncompressed BAM stream and 1--9 for fastest to best compression ratio. Defaults to `4`.

  Please note that both `u` and `0` generates uncompressed output. However, `0` generates BAM stream with `zlib` wrapping while `u` generates raw BAM stream.

Please refer to [`SAMv1.pdf`](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information on SAM/BAM format.

### Headless SAM/BAM Format (`--o-hl_sam`)

This format is specially designed for the `stream` reference parser that allows handling of numerous contigs. As a side effect, the contig name and length information will not be written to SAM header. Each read will be written as unaligned with coordinate information populated in the `OA` tag.

Example:

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt [...]
NM_069135:art_modern:1:nompi:0	4	*	1	0	*	*	1	0	AGCCAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTCGAA	cddddhhhhhhhghhhhhhhghhhhhhhRhhhhhhhghRhhghhh[hhhPhehhhhhhh[hhhfhhhhhdhhhhhhhhfhhe]ehhhhhheg_hhhhQhh[ghhhhhhhhhhdhOhfhhhhhhhh	OA:Z:NM_069135,1,+,3=1D4=1I41=1X27=1I5=1X36=1D5=,255,6;	MD:Z:3^A45T32G36^T5	NM:i:6
[...]
```

This output writer supports other BAM formatting parameters.

## ART-Specific Parameters

- `--id`: Read ID Prefix. Used to distinguish between different runs.
- `--read_len`: Read length. Please note that the read length needs to be smaller than the maximum length supported by the quality profiles.
- `--qual_file_1` and `--qual_file_2`: Quality distribution files for read 1 and read 2 respectively. The first parameter is required while the second parameter is required only for `pe` and `mp` library construction mode.
- `--q_shift_1` and `--q_shift_2`: Shift the base quality of the quality distribution by the specified value. Please note that the shifted quality distribution will be bounded to `--min_qual` and `--max_qual`.
- `--min_qual` and `--max_qual`: Clip the quality distribution to the specified range.
- `--sep_flag`: Separate quality distributions for different bases in the same position.
  - **unset (DEFAULT)**: Use the same quality distribution for different bases in the same position. For example, the quality distribution of pos 10 will be the same regardless whether the incoming base is `A` or `T`.
  - set: Use different quality distributions for different bases in the same position. For example, the quality distribution of pos 10 may **NOT** be the same regardless whether the incoming base is `A` or `T`.
- `--ins_rate_1`, `--ins_rate_2`, `--del_rate_1`, and `--del_rate_2`: Targeted insertion and deletion rate for read 1 and 2.
- `--max_indel`: The maximum total number of insertions and deletions per read.

Quality distributions bundled with the original ART can be found [here](../data/Illumina_profiles).

## Performance Hint

When building `art_modern`, set `USE_HTSLIB` to the latest HTSLib available on your system.  Please also make sure that your HTSLib has been compiled with `-O3 -mtune=native` and linked with [libdeflate](https://github.com/ebiggers/libdeflate). Set `CMAKE_BUILD_TYPE` to `Release` or `RelWithDebInfo`, and `USE_RANDOM_GENERATOR` to `ONEMKL` on Intel/AMD machines.

Also, when executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.

## Generate Illumina Profile

See [ART_profiler_illumina](../deps/ART_profiler_illumina).
