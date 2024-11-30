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
