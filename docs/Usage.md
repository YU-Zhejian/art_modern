# Usage

## Simulation Modes (`--mode`)

## Library Construction Methods (`--lc`)

## Input Formats (`-i-*`)

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

## Output Formats

### Pairwise Alignment Format (`--o-pwa`)

PWA is a serialization of the internally used data structure of the simulator. It is relatively fast and can be converted to other output formats.

### FASTQ Format (`--o-fastq`)

### SAM/BAM Format (`--o-sam`)

### Headless SAM/BAM Format (`--o-hl_sam`)

## ART-Specific Parameters

### Read ID Prefix (`--id`)

### Quality Distribution File (`--qual_file_?`, `--q_shift_?`, `--min_qual`, and `--max_qual`)

### Indel Rate (`--ins_rate_?`, `--del_rate_?`, `--max_indel`)

### Separated Quality Profile for Different Bases (`--sep_flag`)

### Read Length (`--read_len`)

### Paired-End Parameters (`--pe_frag_dist_mean`, `--pe_frag_dist_std_dev`)

## Performance Hint

When building `art_modern`, set `USE_HTSLIB` to the latest HTSLib available on your system.  Please also make sure that your HTSLib has been compiled with `-O3 -mtune=native` and linked with [libdeflate](https://github.com/ebiggers/libdeflate). Set `CMAKE_BUILD_TYPE` to `Release` or `RelWithDebInfo`, and `USE_RANDOM_GENERATOR` to `ONEMKL` on Intel/AMD machines.

Also, when executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.
