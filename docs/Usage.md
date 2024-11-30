# Usage

## Input Formats

Currently, we support input in FASTA and PBSim3 transcripts format.

**FOR FASTA FORMAT**: For read names, only characters before blank space are read.

A compatibility matrix is as follows:

| Parser \ Mode | `wgs`     | `trans`                     | `templ`                     |
|---------------|-----------|-----------------------------|-----------------------------|
| `memory`      | FASTA     | FASTA \| PBSim3 Transcripts | FASTA \| PBSim3 Transcripts |
| `htslib`      | FASTA     | **ERROR**                   | **ERROR**                   |
| `stream`      | **ERROR** | FASTA \| PBSim3 Transcripts | FASTA \| PBSim3 Transcripts |

## Library Construction Methods

## FASTA Parsers

## Performance Hint

Set `USE_HTSLIB` to the latest HTSLib available on your system.  Please also make sure that your HTSLib has been linked with [libdeflate](https://github.com/ebiggers/libdeflate). Set `CMAKE_BUILD_TYPE` to `Release` or `RelWithDebInfo`, and `USE_RANDOM_GENERATOR` to `ONEMKL` on Intel/AMD machines.

Also, when executing `art_modern`, please use `memory` for FASTA parser. Use solid state drive (SSDs) whenever possible. Also use as fewer output writers as possible.
