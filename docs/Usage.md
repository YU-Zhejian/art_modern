# Usage

## Simulation Modes (`--mode`)

## Library Construction Methods (`--lc`)

## Input Formats (`--i-*`)

Currently, we support input in FASTA and PBSIM3 Transcripts format.

**FOR FASTA FORMAT**: For read names, only characters before blank space are read.

A compatibility matrix is as follows:

| Parser \ Mode | `wgs`     | `trans`                     | `template`                  |
|---------------|-----------|-----------------------------|-----------------------------|
| `memory`      | FASTA     | FASTA \| PBSIM3 Transcripts | FASTA \| PBSIM3 Transcripts |
| `htslib`      | FASTA     | **ERROR**                   | **ERROR**                   |
| `stream`      | **ERROR** | FASTA \| PBSIM3 Transcripts | FASTA \| PBSIM3 Transcripts |

`--i-parser`
`--i-type`
`--i-batch_size`
`--i-file`
`--i-fcov`

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

The good old FASTQ format. The qualities are phread encoded in ASCII with an offset of 33.

```text
@NM_069135:art_modern:1:nompi:0
AGCCAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTCGAA
+
BCCCCGGGGGGGFGGGGGGGFGGGGGGG1GGGGGGGFG1GGFGGG:GGG/GDGGGGGGG:GGGEGGGGGCGGGGGGGGEGGD<DGGGGGGDF>GGGG0GG:FGGGGGGGGGGCG.GEGGGGGGGG
```

### SAM/BAM Format (`--o-sam`)

This writes canonical SAM/BAM format using [HTSLib](https://github.com/samtools/htslib). The generated SAM format can be parsed using [samtools](https://github.com/samtools/samtools), and can be used as ground-truth when benchmarking sequence aligners.

This output writer supports computing `NM` and `MD` tag. The mapping qualities for all reads are set to 255.

Example:

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt [...]
@SQ	SN:NM_069135	LN:1718
[...]
NM_069135:art_modern:1:nompi:0	0	NM_069135	1	255	3=1D4=1I41=1X27=1I5=1X36=1D5=	=	1	0	AGCCAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTCGAA	cddddhhhhhhhghhhhhhhghhhhhhhRhhhhhhhghRhhghhh[hhhPhehhhhhhh[hhhfhhhhhdhhhhhhhhfhhe]ehhhhhheg_hhhhQhh[ghhhhhhhhhhdhOhfhhhhhhhh	MD:Z:3^A45T32G36^T5	NM:i:6
[...]
```

This output writer supports BAM format (with default compression rate) and computing CIGAR string using `M` (`BAM_CMATCH`) instead of `=`/`X` (`BAM_CEQUAL`/`BAM_CDIFF`).

### Headless SAM/BAM Format (`--o-hl_sam`)

This format is specially designed for the `stream` reference parser that allows handling of numerous contigs. As a side effect, the contig name and length information will not be written to SAM header. Each read will be written as unaligned with coordinate information populated in the `OA` tag.

Example:

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern --qual_file_1 data/Illumina_profiles/HiSeq2500L125R1.txt [...]
NM_069135:art_modern:1:nompi:0	4	*	1	0	*	*	1	0	AGCCAAACGGGCAACCAGACTCCGCCCATTTCTCAACTCTCTAAGTACCCTGAGAGGTAGTTAGAGAAAACGAGAAACACTACAAACGAGATCAGTTTATCAGCCAATCTGTGTGTTTGTTCGAA	cddddhhhhhhhghhhhhhhghhhhhhhRhhhhhhhghRhhghhh[hhhPhehhhhhhh[hhhfhhhhhdhhhhhhhhfhhe]ehhhhhheg_hhhhQhh[ghhhhhhhhhhdhOhfhhhhhhhh	OA:Z:NM_069135,1,+,3=1D4=1I41=1X27=1I5=1X36=1D5=,255,6;	MD:Z:3^A45T32G36^T5	NM:i:6
[...]
```

This output writer also supports BAM format and computing CIGAR string using `M` (`BAM_CMATCH`) instead of `=`/`X` (`BAM_CEQUAL`/`BAM_CDIFF`).

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
