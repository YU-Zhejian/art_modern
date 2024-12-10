# Usage

## Reading the Documentation

In this documentation, commands will be represented as `ls -lFh` with in-line or block omission represented as `[...]`.

## Simulation Modes (`--mode`)

## Library Construction Methods (`--lc`)

## Input (`--i-*`)

Currently, we support input in FASTA and PBSIM3 Transcripts format. They are controlled by the following major parameters:

- `--i-file`: The input reference file path.
- `--i-type`: The file type of input reference sequences. Currently, FASTA and PBSIM3 Transcripts format are supported.
  - **`auto` (DEFAULT) for extension-based decision.**
    - If the file ends with `.fna`, `.fsa`, `.fa`, `.fasta`, resolve to `fasta`.
    - Otherwise, an error will be raised.
  - `fasta` for FASTA files.
  - `pbsim3_transcripts` for PBSIM3 Transcripts format.

Following are detailed constrains of the aforementioned format:

### Input Parser (`--i-parser`)

- **`auto` (DEFAULT) for size-based determination.**
  - If the file size larger than 1GiB or can not be told (which is quite common if the input was redirected from stdin or other devices), resolve to `htslib` (`wgs` mode) or `stream` (`trans` or `template` mode).
  - Otherwise, use `memory`.
- `memory`: The entire file will be read into the memory.
  - Fast for small reference files.
- `htslib`: Store FASTA Index in memory and fetch sequences from disk using [HTSLib](https://github.com/samtools/htslib).
  - A FASTA Index (Usually ends with `.fai` and can be built using `samtools faidx`) will be needed.
  - Memory-efficient for large genome assembly FASTA files with limited number of long contigs.
  - Inefficient for transcriptome/template FASTAs with large number of (relatively) short contigs.
  - Each thread will hold its own FASTA Index in memory.
  - **DO NOT SUPPORT PBSIM3 TRANSCRIPT FORMAT.**
- `stream`: Streamline the input reference as microbatches and process them one by one.
  - Efficient for transcriptome/template FASTAs with large number of (relatively) short contigs.
  - **CONTIG NUMBER AND LENGTH INFORMAION NOT AVAILABLE**, so use headless SAM/BAM output if a SAM/BAM output is needed.
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

### Conclusive Remarks

A compatibility matrix is as follows:

| Parser \ Mode | `wgs`     | `trans`                     | `template`                  |
|---------------|-----------|-----------------------------|-----------------------------|
| `memory`      | FASTA     | FASTA \| PBSIM3 Transcripts | FASTA \| PBSIM3 Transcripts |
| `htslib`      | FASTA     | **ERROR**                   | **ERROR**                   |
| `stream`      | **ERROR** | FASTA \| PBSIM3 Transcripts | FASTA \| PBSIM3 Transcripts |

## Output Formats (`--o-*`)

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

### Read ID Prefix (`--id`)

The prefix of read IDs. Used to distinguish between different runs.

### Quality Distribution File (`--qual_file_?`, `--q_shift_?`, `--min_qual`, and `--max_qual`)

### Indel Rate (`--ins_rate_?`, `--del_rate_?`, `--max_indel`)

### Separated Quality Profile for Different Bases (`--sep_flag`)

### Read Length (`--read_len`)

### Paired-End Parameters (`--pe_frag_dist_mean`, `--pe_frag_dist_std_dev`)

## Performance Hint

When building `art_modern`, set `USE_HTSLIB` to the latest HTSLib available on your system.  Please also make sure that your HTSLib has been compiled with `-O3 -mtune=native` and linked with [libdeflate](https://github.com/ebiggers/libdeflate). Set `CMAKE_BUILD_TYPE` to `Release` or `RelWithDebInfo`, and `USE_RANDOM_GENERATOR` to `ONEMKL` on Intel/AMD machines.

Also, when executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.
