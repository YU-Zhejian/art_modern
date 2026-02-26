# Usage of `art_modern`

This is the detailed usage of `art_modern`. Every parameter and their combinations are introduced in detail. In this documentation, commands will be represented as `ls -lFh` with in-line or block omission represented as `[...]`.

## Parameter Controlling Simulation Behavior

### Simulation Modes (`--mode`)

The simulator currently supports 3 simulation modes:

- **`wgs` (DEFAULT) for whole-genome sequencing.** For scenarios with constant coverage over a limited number of sized contigs.
- `trans` for transcriptome simulation. For scenarios with many short contigs and/or contig-specific coverage.
- `template` for template simulation. Often being used with a high-level simulator that produces read fragments.

### Library Construction Methods (`--lc`)

The simulator currently supports 3 library construction methods:

#### **`se` (DEFAULT) for single-end reads**

- **For `wgs` and `trans` mode:** A read may start from any position in the reference contig that is long enough to contain the read.
- **For `template` mode:** A read will start from position 0 in the reference contig.

A diagram for single-end read extraction:

```text
FRAGMENT:  |==============================|
SE READ:   |------->
```

#### `pe` for paired-end reads

- **For `wgs` and `trans` mode:** A fragment (called insert in some other simulators), whose length is drawn from a Gaussian distribution specified by `--pe_frag_dist_mean` and `--pe_frag_dist_std_dev`, is created from any position in the reference contig that is long enough to contain the fragment.
- **For `template` mode:** A fragment that spans the entire length of the reference contig is created.
- Paired-end reads facing inward of the fragment are then created from both ends of the fragment.

A diagram for paired-end read extraction:

```text
FRAGMENT:  |==============================|
PE READ 1: |------->
PE READ 2:                        |<------|
```

#### `mp` for mate-paired reads

- Extraction of the fragment as-is in `pe` mode.
- Mate-paired reads facing outward of the fragment are created from both ends of the fragment.

A diagram for mate-pair extraction:

```text
FRAGMENT:  |==============================|
MP READ 1:                        |------->
MP READ 2: <-------|
```

#### See also

- [Illumina paired vs. single end blog](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html).
- [Illumina explaination of mate-pair simulation](https://www.illumina.com/science/technology/next-generation-sequencing/mate-pair-sequencing.html).

(parallelism-section)=
### Thread-Level Parallelism (`--parallel`)

Number of threads to use. Use `-1` to disable parallelism. Use `0` (default) to use all available cores. Oversubscription is allowed but not recommended.

**NOTE** If MPI is enabled, each MPI process will use the specified number of threads. For example, if you run `mpiexec -n 4 art_modern --parallel 2 ...`, there will be 8 threads in total.

**NOTE** SAM/BAM output writer is parallelized using its own thread pool. See parameter [`--o-sam-num_threads`](#out-sam-section) and [`--o-hl_sam-num_threads`](#out-hl_sam-section) for more information. For example, if you run `art_modern --parallel 4 --o-sam-num_threads 2 --o-hl_sam-num_threads 3 ...`, there will be 9 threads in total.

**WARNING** Do not oversubscribe the number of threads. It impacts the performance of the program.

### Error-Profile Related Parameters

#### Which Quality Distribution File to Use

Quality distribution files specify the base quality distribution at each position of the read. They are used to simulate sequencing errors.

Users may supply self-trained quality distribution files using `--qual_file_1` and `--qual_file_2` for read 1 and read 2 respectively. The format of quality distribution files is compatible with the original ART. The first parameter is required while the second parameter is required only for `pe` and `mp` library construction mode.

Users may also use built-in quality distribution files using `--builtin_qual_file` (added in [1.1.0](#v-1.1.0-section).). The built-in quality distribution files are bundled with the simulator. Note that the built-in quality distribution files comes with different maximum read-length and different support over paired-end reads. So please make sure that the read length specified by `--read_len` is compatible with the built-in quality distribution file.

Quality distributions bundled with the original ART can be found [here](https://github.com/YU-Zhejian/art_modern/tree/master/data/Illumina_profiles). In details:

Here lists existing Illumina models supported by the original ART, which is still available in `art_modern`.

| Name   | Name (`art_modern`)     | R1 `txt` file name    | R2? | RLen1 | RLen2 | Model Name        |
|--------|-------------------------|-----------------------|-----|-------|-------|-------------------|
| `GA1`  | `GA1Recalibrated_36bp`  | `EmpR36R1`            | Y   | 36    | 36    | GenomeAnalyzer I  |
| `GA1`  | `GA1Recalibrated_44bp`  | `EmpR44R1`            | Y   | 44    | 44    | GenomeAnalyzer I  |
| `GA2`  | `GA2Recalibrated_50bp`  | `EmpR50R1`            | Y   | 50    | 50    | GenomeAnalyzer II |
| `GA2`  | `GA2Recalibrated_75bp`  | `EmpR75R1`            | Y   | 75    | 75    | GenomeAnalyzer II |
| `HS10` | `HiSeq1000_100bp`       | `Emp100R1`            | Y   | 100   | 100   | HiSeq 1000        |
| `MSv1` | `MiSeq_250bp`           | `EmpMiSeq250R1`       | Y   | 250   | 250   | MiSeq             |
| `HS20` | `HiSeq2000_100bp`       | `HiSeq2kL100R1`       | Y   | 100   | 100   | HiSeq 2000        |
| `HS25` | `HiSeq2500_125bp`       | `HiSeq2500L125R1`     | Y   | 126   | 126   | HiSeq 2500        |
| `HS25` | `HiSeq2500_150bp`       | `HiSeq2500L150R1`     | Y   | 150   | 150   | HiSeq 2500        |
| `HSXn` | `HiSeqX_PCR_Free_150bp` | `HiSeqXPCRfreeL150R1` | Y   | 151   | 151   | HiSeqX PCR free   |
| `HSXt` | `HiSeqX_TruSeq_150bp`   | `HiSeqXtruSeqL150R1`  | Y   | 151   | 151   | HiSeqX TruSeq     |
| `MinS` | `MiniSeq_TruSeq_50bp`   | `MiniSeqTruSeqL50`    | N   | 51    | N/A   | MiniSeq TruSeq    |
| `MSv3` | `MiSeq_v3_250bp`        | `MiSeqv3L250R1`       | Y   | 251   | 251   | MiSeq v3          |
| `NS50` | `NextSeq500_v2_75bp`    | `NextSeq500v2L75R1`   | Y   | 76    | 76    | NextSeq500 v2     |

Following profiles are only available in `art_modern` but not in the original ART as a built-in model. However, the
model error files are provided in the original ART.

| Name (`art_modern`)       | R1 `txt` file name      | R2? | RLen1 | RLen2 | Model Name        |
|---------------------------|-------------------------|-----|-------|-------|-------------------|
| `GA1_36bp`                | `Emp36R1`               | Y   | 36    | 36    | GenomeAnalyzer I  |
| `GA1_44bp`                | `Emp44R1`               | Y   | 44    | 44    | GenomeAnalyzer I  |
| `GA2_50bp`                | `Emp50R1`               | Y   | 53    | 50    | GenomeAnalyzer II |
| `GA2_75bp`                | `Emp75R1`               | Y   | 75    | 75    | GenomeAnalyzer II |
| `HiSeq2500Filtered_150bp` | `HiSeq2500L150R1filter` | Y   | 150   | 150   | HiSeq 2500        |

#### Shift and Clip Parameters

`--q_shift_1` and `--q_shift_2`: Shift the base quality of the quality distribution by the specified value. Please note that the shifted quality distribution will be bounded to `--min_qual` and `--max_qual`.

`--min_qual` and `--max_qual`: Clip the quality distribution to the specified range.

#### Separate Quality Distributions for Different Bases (`--sep_flag`)

This flag can separate quality distributions for different bases in the same position.

- **unset (DEFAULT)**: Use the same quality distribution for different bases in the same position. For example, the quality distribution of pos 10 will be the same regardless whether the incoming base is `A` or `T`.
- set: Use different quality distributions for different bases in the same position. For example, the quality distribution of pos 10 may **NOT** be the same regardless whether the incoming base is `A` or `T`.

## Parameters Controlling Input

### Input File Path (`--i-file`)

The input reference file path. It must exist in the filesystem. Device paths like `/dev/stdin` is supported for non-MPI builds.

Currently, we support input in FASTA and PBSIM3 Transcripts format. An additional parameter [`--i-type`](#input-file-type-section) is used to specify the file type.

**NOTE** Some parser may support compressed FASTAs. See [`--i-parser`](#input-parser-section) for more detils.

(input-file-type-section)=
### Input File Type (`--i-type`)

The file type of input reference sequences. Currently, FASTA and PBSIM3 Transcripts format are supported.

- **`auto` (DEFAULT) for extension-based decision.**
  - If the file ends with `.fna`, `.fsa`, `.fa`, `.fasta`, resolve to `fasta`.
  - Otherwise, an error will be raised.
- `fasta` for FASTA files.
- `pbsim3_transcripts` for PBSIM3 Transcripts format.

**NOTE** If you're using UNIX devices like `/dev/stdin` as input, please make sure that this option is specified.

(input-parser-section)=
### Input Parser (`--i-parser`)

#### **`auto` (DEFAULT) for size-based determination.**

If the file size larger than 1 GiB or cannot be told (which is quite common if the input was redirected from stdin or other devices), resolve to `htslib` (`wgs` mode) or `stream` (`trans` or `template` mode). Otherwise, use `memory`.

#### `memory` for in-memory reference parser

The entire file will be read into the memory. Fast for small reference files.

#### `htslib` for HTSLib FAIDX parser

Store FASTA Index in memory and fetch sequences from disk using [HTSLib](https://github.com/samtools/htslib). Memory-efficient for large genome assembly FASTA files with limited number of long contigs, but inefficient for transcriptome/template FASTAs with large number of (relatively) short contigs.

The FASTA index, usually ends with `.fai` and can be built using `samtools faidx`, will be created automatically if not exist.

**NOTE** For `wgs` mode only.

**NOTE** Do not support PBSIM3 Transcripts format.

**NOTE** Do not support UNIX devices for input.

**NOTE** `htslib` parser can be used on a `.fa.gz` file. However, please note that:

- The GZip file **MUST** be [`bgzip`](https://www.htslib.org/doc/bgzip.html)-compressed. Plain GZip will **NOT** work. If you try to use a GZip-compressed FASTA file, you'll see error messages like:

  ```text
  [E::fai_build_core] File truncated at line 1
  [E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
  [E::fai_path] Failed to build index file for reference file [...]
  ```

- Specify `--i-type fasta` to use this parser.

#### `stream` for one-pass streaming parser

Streamline the input reference as batches and process them one by one. Efficient for transcriptome/template FASTAs with large number of (relatively) short contigs.

Please note that here, the contig number and length are not accessible to the output dispatcher. So use [headless SAM/BAM output](#out-hl_sam-section) if a SAM/BAM output is needed.

The batch size of each batch is controlled by `--i-batch_size`. Larger batch size may improve performance but requires more memory. Use larger batch size if the targeted sequencing depth is relatively shallow.

### Coverage (`--i-fcov`)

This option allows you to specify coverage. `art_modern` supports the following coverage mode:

- Unified coverage floating point in `double` data type (e.g., `10.0`). Under this scenario, the coverage is identical for all contigs.
- Per-contig coverage headless tab-separated value (TSV) without strand information.
- Per-contig coverage headless TSV with strand information.

**NOTE** This parameter is ignored if the file type is set to `pbsim3_transcripts`. This format comes with contig-specific coverage information.

**NOTE** We prefer unified coverage floating points to TSVs. That is, if you set this parameter to 10.0, we're **NOT** going to check whether there's a file named `10.0`. So please make sure that the file name is not a number if you want to specify a coverage file.

**NOTE** For coverage of the template mode:

- If a unified coverage file is provided, the coverage will be intepreted as positive coverage. That is:
  - If the library is SE, all reads will be generated from the positive strand at the 1st base.
  - If the library is PE, all read 1 will be generated from the positive strand at the first base, with read 2 generated from the negative strand at the last base.
  - If the library is MP, all read 1 will be generated from the positive strand to the last base, with read 2 generated from the negative strand to the first base.
- If a 2-column (unstranded) coverage file is provided, the coverage will be intepreted as positive coverage.
  - Introduced in [1.3.0](#v-1.3.0-section).
- If a 3-column (stranded) coverage file or input in format of `pvsim3_transcripts` is provided, the coverage will be intepreted as-is.

For example, if you specified `--read_len_1 10 --read_len_2 150` with unified/unstranded coverage 5 and PE library construction mode, you'll be likely to see:

```text
|================================================|
|---->                              <------------|
|---->                              <------------|
|---->                              <------------|
|---->                              <------------|
|---->                              <------------|
```

With stranded coverage positive=3 and negative=2, you may see:

```text
|================================================|
|---->                              <------------|
|---->                              <------------|
|---->                              <------------|
|------------>                              <----|
|------------>                              <----|
```

while for MP library construction mode you would see:

```text
|================================================|
<------------|                              |---->
<------------|                              |---->
<------------|                              |---->
<----|                              |------------>
<----|                              |------------>
```

### Compatibility Matrices of Input Parameters

Compatibility matrix of file type, simulation mode, and parser:

| Parser \ Mode | `wgs`     | `trans`                        | `template`                     |
|---------------|-----------|--------------------------------|--------------------------------|
| `memory`      | `fasta`   | `fasta` / `pbsim3_transcripts` | `fasta` / `pbsim3_transcripts` |
| `htslib`      | `fasta`   | **ERROR**                      | **ERROR**                      |
| `stream`      | **ERROR** | `fasta` / `pbsim3_transcripts` | `fasta` / `pbsim3_transcripts` |

Compatibility matrix of coverage mode, simulation mode, and file type:

| Mode \ File Type | `fasta`                        | `pbsim3_transcripts` |
|------------------|--------------------------------|----------------------|
| `wgs`            | Unified                        | **ERROR**            |
| `trans`          | Unified / Unstraded / Stranded | **IGNORED**          |
| `template`       | Unified / Unstraded / Stranded | **IGNORED**          |

## Parameters Controlling Output

Here introduces diverse output formats supported by `art_modern`. You may specify none of them to perform simulation without any output for benchmarking purposes.

If the parent directory does not exist, it will be created using a `mkdir -p`-like manner.

For all output format, a queue size parameter (e.g., `--o-pwa-queue_size`) is provided to control the size of lockfree queue. Turn this up to improve performance when you have enough memory, and low if the application is consuming too much memory. Added in [1.3.2](#v-1.3.2-section).

### Pairwise Alignment Format (`--o-pwa`)

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

This output format supports [`am_compress`](#am_compress-usage-section)-like compression arguments (Namely, `--o-pwa-compression`, `--o-pwa-compression_level`, `--o-pwa-buffer_size` and `--o-pwa-num_threads`). Added in [1.3.4](#v-1.3.4-section).

### FASTQ Format (`--o-fastq`)

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

For PE/MP simulation: The generated file is interleaved FASTQ. Split to 2 FASTQ files using [`seqtk`](https://github.com/lh3/seqtk):

```shell
# Read 1
seqtk seq a.fq -1 > a_1.fq
# Read 2
seqtk seq a.fq -2 > a_2.fq
```

This output format supports [`am_compress`](#am_compress-usage-section)-like compression arguments. Added in [1.3.4](#v-1.3.4-section).

See also:

- [Specifications of Common File Formats Used by the ENCODE Consortium at UCSC](https://genome.ucsc.edu/ENCODE/fileFormats.html#FASTQ).
- [Common File Formats Used by the ENCODE Consortium](https://www.encodeproject.org/help/file-formats/#fastq).

### FASTA Format (`--o-fasta`)

Sequence-only output format. Example:

```text
>NM_069135:art_modern:1:nompi:0
AGCCAAACGGGCAACCAGACTCCGCC[...]
```

This simulator generates one-line FASTA files. That is, each sequence occupies only one line.

This output format supports [`am_compress`](#am_compress-usage-section)-like compression arguments. Added in [1.3.4](#v-1.3.4-section).

(out-sam-section)=
### SAM/BAM Format (`--o-sam`)

Sequence Alignment/Map (SAM) and Binary Alignment/Map (BAM) format supports storing of ground-truth alignment information and other miscellaneous parameters. They can be parsed using [samtools](https://github.com/samtools/samtools), [`pysam`](https://pysam.readthedocs.io/) and other libraries, and can be used as ground-truth when benchmarking sequence aligners.

This writes canonical SAM/BAM format using [HTSLib](https://github.com/samtools/htslib). This output writer supports computing `NM` and `MD` tag. The mapping qualities for all reads are set to 255.

Example (Tabs replaced with spaces):

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern [...]
@SQ	SN:NM_069135	LN:1718
[...]
NM_069135:[...] 0 NM_069135 1 255 3=[...]5= = 1 0 AGCC[...] BCCC[...] MD:[...] NM:[...]
[...]
```

Other SAM/BAM formatting parameters includes:

- `--o-sam-use_m`: Computing CIGAR string using `M` (`BAM_CMATCH`) instead of `=`/`X` (`BAM_CEQUAL`/`BAM_CDIFF`). Rarely used in next-generation sequencing but common for long-read sequencing alignments.
- `--o-sam-write_bam`: Write BAM instead of SAM.
- `--o-sam-num_threads`: Number of threads used to compress BAM output.
- `--o-sam-compress_level`: [`zlib`](https://www.zlib.net/) compression level. Supports `[u0-9]` with `u` for uncompressed BAM stream, 0 for uncompressed stream with `zlib` wrapping, and 1--9 for fastest to best compression ratio. Defaults to `4`. Please note that both `u` and `0` generates uncompressed output. However, `0` generates BAM stream with `zlib` wrapping while `u` generates raw BAM stream.
- `--o-sam-without_tag_MD`: Do not write `MD` tag. Added in [1.3.3](#v-1.3.3-section).
- `--o-sam-without_tag_NM`: Do not write `NM` tag. Added in [1.3.3](#v-1.3.3-section).
- `--o-sam-no_qual`: Do not write quality string (Write quality string as `*`). Added in [1.3.3](#v-1.3.3-section).

**NOTE** For PE/MP simulation, the 2 reads from one fragment may **NOT** appear together.

See also: [`SAMv1.pdf`](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information on SAM/BAM format.

(out-hl_sam-section)=
### Headless SAM/BAM Format (`--o-hl_sam`)

This format is specially designed for the `stream` reference parser that allows handling of numerous contigs. As a side effect, the contig name and length information will not be written to SAM header. Each read will be written as unaligned with coordinate information populated in the `OA` tag.

Example (Tabs replaced with spaces):

```text
@HD	VN:1.4	SO:unsorted
@PG	ID:01	PN:art_modern	CL:opt/build_debug/art_modern [...]
NM_069135:[...] 4 * 1 0 * * 1 0 AGCC[...] BCCC[...] OA:[...] MD:[...] NM:[...]
[...]
```

An additional option: `--o-hl_sam-without_tag_OA`, setting such which will not write `OA` tag. Added in [1.3.3](#v-1.3.3-section).

This output writer supports other BAM formatting parameters.

## Miscellaneous Parameters

### Read Length Parameters

`--read_len`, `--read_len_1`, and `--read_len_2` controls read length. Specifically:

- If none of the 3 parameters are set, will use the longest read length supported by the quality distribution file(s).
- If only `--read_len` is set, will use the specified read length for both read 1 and read 2.
- `--read_len` cannot be used together with `--read_len_1` or `--read_len_2`.
- If `--read_len_1` and/or `--read_len_2` is set, will use the specified read lengths for read 1 and/or read 2 respectively.
- For paired-end or mate-pair library construction mode, both read lengths must be set either using `--read_len` **OR** `--read_len_1` and `--read_len_2`.

### Indel Error Parameters

- `--ins_rate_1`, `--ins_rate_2`, `--del_rate_1`, and `--del_rate_2`: Targeted insertion and deletion rate for read 1 and 2.
- `--max_indel`: The maximum total number of insertions and deletions per read.

### Ambiguous Base Parameters

- `--max_n`: The maximum total number of ambiguous bases (N) per read.

### Others

- `--id`: Read ID Prefix. Used to distinguish between different runs.
- `--i-seed`: Seed for random number generator to ensure reproducibility. Added in [1.3.2](#v-1.3.2-section). **NOTE** To ensure reproducibility, please use the same number of threads, number of MPI progress under the same random generator when running `art_modern` with the same seed.

## Logging Parameters

- `--reporting_interval-job_executor`: The reporting interval (in seconds) for individual `art_modern` job. Added in [1.2.1](#v-1.2.1-section).
- `--reporting_interval-job_pool`: The reporting interval (in seconds) for ThreadPool. Added in [1.2.1](#v-1.2.1-section).

Those logs are designed for observation of long-running jobs. Increase the interval to reduce log size.

(am-environment-variables-section)=
## Environment Variables

- `ART_NO_LOG_DIR`: If set (to any value), disables creation of a log directory. Logging to files is skipped, and a warning is printed. Added in [1.1.4](#v-1.1.4-section).
  - **NOTE** If you're using `art_modern` with MPI, this discards all logs generated from rank 1, 2, 3....
- `ART_LOG_DIR`: If set, specifies the directory where log files will be written. If not set, defaults to `log.d` and a warning is printed. Ignored if `ART_NO_LOG_DIR` is set. Added in [1.1.4](#v-1.1.4-section).
  - **NOTE** If you're running 2 `art_modern` processes simultaneously, please make sure that they are using different log directories.
