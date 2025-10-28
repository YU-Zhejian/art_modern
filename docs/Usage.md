# Usage

This is the detailed usage of `art_modern` and its companying binaries. Every parameter and their combinations are introduced in detail. In this documentation, commands will be represented as `ls -lFh` with in-line or block omission represented as `[...]`.

## The Main Simulator: `art_modern` Executable

### Simulation Modes (`--mode`)

- **`wgs` (DEFAULT) for whole-genome sequencing.** For scenarios with constant coverage over a limited number of sized contigs.
- `trans` for transcriptome simulation. For scenarios with many short contigs and/or contig-specific coverage.
- `template` for template simulation. Often being used with a high-level simulator that produces read fragments.

### Library Construction Methods (`--lc`)

See [here](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) and [here](https://www.illumina.com/science/technology/next-generation-sequencing/mate-pair-sequencing.html) for details on library construction methods.

- **`se` (DEFAULT) for single-end reads.** That includes:
  - **For `wgs` and `trans` mode:** A read may start from any position in the reference contig that is long enough to contain the read.
  - **For `template` mode:** A read will start from position 0 in the reference contig.
- `pe` for paired-end reads. That includes:
  - **For `wgs` and `trans` mode:** A fragment, whose length is drawn from a Gaussian distribution specified by `--pe_frag_dist_mean` and `--pe_frag_dist_std_dev`, is created from any position in the reference contig that is long enough to contain the fragment.
  - **For `template` mode:** A fragment that spans the entire length of the reference contig is created.
  - Paired-end reads facing inward of the fragment are then created from both ends of the fragment.
- `mp` for mate-paired reads. That includes:
  - Extraction of the fragment as-is in `pe` mode.
  - Mate-paired reads facing outward of the fragment are created from both ends of the fragment.

Parameters `--pe_frag_dist_mean` and `--pe_frag_dist_std_dev` are needed to specify the length distribution of fragments. In original ART, the fragment lengths obey Gaussian (normal) distribution.

(parallelism-section)=
### Parallelism (`--parallel`)

Number of threads to use. Use `-1` to disable parallelism. Use `0` (default) to use all available cores.

### Input (`--i-*`)

Currently, we support input in FASTA and PBSIM3 Transcripts format. They are controlled by the following major parameters:

- `--i-file`: The input reference file path. It must exist in the filesystem. Device paths like `/dev/stdin` is supported.
- `--i-type`: The file type of input reference sequences. Currently, FASTA and PBSIM3 Transcripts format are supported.
  - **`auto` (DEFAULT) for extension-based decision.**
    - If the file ends with `.fna`, `.fsa`, `.fa`, `.fasta`, resolve to `fasta`.
    - Otherwise, an error will be raised.
  - `fasta` for FASTA files.
  - `pbsim3_transcripts` for PBSIM3 Transcripts format.

#### Input Parser (`--i-parser`)

- **`auto` (DEFAULT) for size-based determination.**
  - If the file size larger than 1 GiB or cannot be told (which is quite common if the input was redirected from stdin or other devices), resolve to `htslib` (`wgs` mode) or `stream` (`trans` or `template` mode).
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

#### Coverage (`--i-fcov`)

**NOTE** This parameter is ignored if the file type is set to `pbsim3_transcripts`.

This option allows you to specify coverage. `art_modern` supports the following coverage mode:

- Unified coverage in `double` data type (e.g., 10.0). Under this scenario, the coverage is identical for all contigs.
- Per-contig coverage without strand information.
- Per-contig coverage with strand information.

**NOTE** We prefer unified coverage. That is, if you set this parameter to 10.0, we're **NOT** going to check whether there's a file named `10.0`. So please make sure that the file name is not a number if you want to specify per-contig coverage.

#### Conclusive Remarks

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

### Output Formats (`--o-*`)

Here introduces diverse output formats supported by `art_modern`. You may specify none of them to perform simulation without any output for benchmarking purposes.

If the parent directory does not exist, it will be created using a `mkdir -p`-like manner.

#### Pairwise Alignment Format (`--o-pwa`)

PWA is serialization of the internally used data structure of the simulator. Currently, PWA consists of some metadata lines (Started with `#`) and a list of alignments, each spanning 4 lines.

The alignment lines are as follows:

1. `>`, read name, a blank, and the leftmost position of the alignment.
2. The aligned read sequence, gapped with `-`.
3. The aligned reference sequence, gapped with `-`.
4. The gapless quality string.

Example:

```text
#PWA
#ARGS: opt/build_debug/art_modern [...]
>NM_069135:art_modern:1:nompi:0	NM_069135:0:+
AGC-CAAACGGGCAACCAGACTCCGC[...]
AGCACAAA-GGGCAACCAGACTCCGC[...]
BCCCCGGGGGGGFGGGGGGGFGGGGG[...]
[...]
```

#### FASTQ Format (`--o-fastq`)

The good old FASTQ format. The qualities are Phread encoded in ASCII with an offset of 33.

```text
@NM_069135:art_modern:1:nompi:0
AGCCAAACGGGCAACCAGACTCCGCC[...]
+
BCCCCGGGGGGGFGGGGGGGFGGGGG[...]
```

FASTQ files can be easily converted to FASTA using [`seqtk`](https://github.com/lh3/seqtk):

```shell
cat in.fq | seqtk seq -A > out.fa
```

See also:

- [Specifications of Common File Formats Used by the ENCODE Consortium at UCSC](https://genome.ucsc.edu/ENCODE/fileFormats.html#FASTQ).
- [Common File Formats Used by the ENCODE Consortium](https://www.encodeproject.org/help/file-formats/#fastq)

#### FASTA Format (`--o-fasta`)

FASTA format output. Note that qualities are not stored in this format. Example:

```text
>NM_069135:art_modern:1:nompi:0
AGCCAAACGGGCAACCAGACTCCGCC[...]
```

#### SAM/BAM Format (`--o-sam`)

Sequence Alignment/Map (SAM) and Binary Alignment/Map (BAM) format supports storing of ground-truth alignment information and other miscellaneous parameters. They can be parsed using [samtools](https://github.com/samtools/samtools), [`pysam`](https://pysam.readthedocs.io/) and other libraries, and can be used as ground-truth when benchmarking sequence aligners.

This writes canonical SAM/BAM format using [HTSLib](https://github.com/samtools/htslib). This output writer supports computing `NM` and `MD` tag. The mapping qualities for all reads are set to 255.

Example:

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern [...]
@SQ	SN:NM_069135	LN:1718
[...]
NM_069135:art_modern:1:nompi:0	0	NM_069135	1	255	3=1D4=1I41=1X27=1I5=1X36=1D5=	=	1	0	AGCCAAACGGGCAACCAGACTCCGCC[...]	BCCCCGGGGGGGFGGGGGGGFGGGGG[...]	MD:Z:3^A45T32G36^T5	NM:i:6
[...]
```

Other SAM/BAM formatting parameters includes:

- `--o-sam-use_m`: Computing CIGAR string using `M` (`BAM_CMATCH`) instead of `=`/`X` (`BAM_CEQUAL`/`BAM_CDIFF`). Rarely used in next-generation sequencing but common for long-read sequencing alignments.
- `--o-sam-write_bam`: Write BAM instead of SAM.
- `--o-sam-num_threads`: Number of threads used to compress BAM output.
- `--o-sam-compress_level`: [`zlib`](https://www.zlib.net/) compression level. Supports `[u0-9]` with `u` for uncompressed BAM stream and 1--9 for fastest to best compression ratio. Defaults to `4`.

  Please note that both `u` and `0` generates uncompressed output. However, `0` generates BAM stream with `zlib` wrapping while `u` generates raw BAM stream.

See also: [`SAMv1.pdf`](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information on SAM/BAM format.

#### Headless SAM/BAM Format (`--o-hl_sam`)

This format is specially designed for the `stream` reference parser that allows handling of numerous contigs. As a side effect, the contig name and length information will not be written to SAM header. Each read will be written as unaligned with coordinate information populated in the `OA` tag.

Example:

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern [...]
NM_069135:art_modern:1:nompi:0	4	*	1	0	*	*	1	0	AGCCAAACGGGCAACCAGACTCCGCC[...]	BCCCCGGGGGGGFGGGGGGGFGGGGG[...]	OA:Z:NM_069135,1,+,3=1D4=1I41=1X27=1I5=1X36=1D5=,255,6;	MD:Z:3^A45T32G36^T5	NM:i:6
[...]
```

This output writer supports other BAM formatting parameters.

### ART-Specific Parameters

- `--id`: Read ID Prefix. Used to distinguish between different runs.
- `--read_len`: Read length. Please note that the read length needs to be smaller than the maximum length supported by the quality profiles.
- `--qual_file_1` and `--qual_file_2`: Quality distribution files for read 1 and read 2 respectively. The first parameter is required while the second parameter is required only for `pe` and `mp` library construction mode. `--builtin_qual_file` may be used to specify a built-in quality distribution.
- `--q_shift_1` and `--q_shift_2`: Shift the base quality of the quality distribution by the specified value. Please note that the shifted quality distribution will be bounded to `--min_qual` and `--max_qual`.
- `--min_qual` and `--max_qual`: Clip the quality distribution to the specified range.
- `--sep_flag`: Separate quality distributions for different bases in the same position.
  - **unset (DEFAULT)**: Use the same quality distribution for different bases in the same position. For example, the quality distribution of pos 10 will be the same regardless whether the incoming base is `A` or `T`.
  - set: Use different quality distributions for different bases in the same position. For example, the quality distribution of pos 10 may **NOT** be the same regardless whether the incoming base is `A` or `T`.
- `--ins_rate_1`, `--ins_rate_2`, `--del_rate_1`, and `--del_rate_2`: Targeted insertion and deletion rate for read 1 and 2.
- `--max_indel`: The maximum total number of insertions and deletions per read.

Quality distributions bundled with the original ART can be found [here](../data/Illumina_profiles).

### Environment Variables

- `ART_NO_LOG_DIR`: If set (to any value), disables creation of a log directory. Logging to files is skipped, and a warning is printed.
  - **NOTE** If you're using `art_modern` with MPI, this discards all logs generated from rank 1, 2, 3....
- `ART_LOG_DIR`: If set, specifies the directory where log files will be written. If not set, defaults to `log.d` and a warning is printed. Ignored if `ART_NO_LOG_DIR` is set.
  - **NOTE** If you're running 2 `art_modern` processes simultaneously, please make sure that they are using different log directories.

## Building New ART Profiles: `art_profile_builder` Executable

This new executable is designed to supress the old `art_profiler_illumina` Shell/Perl scripts for building ART-compatible quality profiles. It supports building profiles from FASTQ and SAM/BAM files.

**NOTE** All reads with quality will be used to build the profile. So if you use SAM/BAM files as input, please make sure that primiary- and un-aligned reads have quality information, and secondary- or supplementary-aligned reads have no quality information.

### Options

- `--i-file`: Input file path. FASTQ/SAM/BAM files are supported. We currently relies on HTSLib to determine file format automatically. See also: [`htsfile(1)`](https://www.htslib.org/doc/htsfile.html).
- `--i-num_threads`: Number of threads used to decompress BAM input.
- `--read_len`: Read length. For each read, only the first `read_len` bases are used. If the read is shorter than `read_len`, only the available bases are used.
- `--is_pe`: For SAM/BAM input, whether the input is paired-end.
- `--o-file1`: Output quality profile file path for read 1.
- `--o-file2`: Output quality profile file path for read 2. Required if `--is_pe` is set.
- `--parallel`: Same as [`art_modern`'s `--parallel`](#parallelism-section).
- `--old_behavior`:  Simulate the behaviour of original ART profile builder. If set, all qualities will be offsetted by 1. See also: [The Bug in `Srcurity.md`](#original-art-profile-builder-bug).

### Example

Building profiles from single-end FASTQ files:

```shell
art_profile_builder \
    --i-file out_se.fq \
    --read_len 36 \
    --o-file1 out_se_art_cxx_fq.txt \
    --parallel 2 \
    --i-num_threads 4
```

Building profiles from paired-end FASTQ files:

```shell
for i in 1 2; do
    art_profile_builder \
        --i-file out_pe."${i}".fq \
        --read_len 36 \
        --o-file1 out_pe_art_cxx_fq_R"${i}".txt \
        --parallel 2 \
        --i-num_threads 4
done
```

Building profiles from single-end SAM/BAM files:

```shell
art_profile_builder \
    --i-file out_se.sam \
    --read_len 36 \
    --o-file1 out_se_art_cxx_sam.txt \
    --parallel 2 \
    --i-num_threads 4
```

Building profiles from paired-end SAM/BAM files:

```shell
art_profile_builder \
    --i-file out_pe.sam \
    --read_len 36 \
    --is_pe \
    --o-file1 out_pe_art_cxx_sam_R1.txt \
    --o-file2 out_pe_art_cxx_sam_R2.txt \
    --parallel 2 \
    --i-num_threads 4
```

## Appendix

### More Instructions on FASTA Format

FASTA format can be parsed by all parsers. However, please keep in mind that for `memory` parser, the following kinds of FASTA files are supported:

```text
>one_line_fasta
AAAAAAAAAAAAAAAAA
>multi_line_fasta
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
>fasta_with_empty_sequence
>fasta_with_empty_sequence_with_newlines



>fasta_with_spaces_in_name some description here
AAAAAAAAA
```

Note that `fasta_with_empty_sequence_with_newlines` is **NOT** supported by [PacBio Formats](https://pacbiofileformats.readthedocs.io/en/13.0/FASTA.html) or [NCBI GenBank FASTA Specification](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) and [NCBI GenBank Submission Guidelines](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#files).

The following kinds of FASTA files are **NOT** supported:

```text
> some_sequence_without_a_name
AAAAAAAA
```

Also, note that all characters other than `ACGTacgt` will be regarded as `N`. We do **NOT** support IUPAC codes.

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

Using empty file or `/dev/null` as input is allowed since 1.1.9. It will generate empty FASTA/FASTQ/PWA files as output. Remember to specify `memory` or `stream` as `--i-parser` `--i_type` and do **NOT** use SAM/BAM output writer in this case (as SAM/BAM output writer will think you're using streamed input and raise an exception).

### Run-Time Performance Hint

When executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.

SAM/BAM output writers are memory- and time-consuming due to compression. If you don't need SAM/BAM output, please don't enable it.
