# Readme

Here records scripts and notes about additional models for `art_modern`.

## Adding New Models

Following is a table of sequencing instrument models that are manufactured by Illumina, BGI/MGI, and PacBio (Onso)
extracted from the EBI ENA database as of December 2025, using `fetch.sh`.

| Instrument Model             | Count      | Supported? |
|------------------------------|------------|------------|
| illumina novaseq 6000        | 10,562,413 | N          |
| illumina miseq               | 9,004,342  | Y          |
| illumina hiseq 2500          | 4,181,585  | Y          |
| illumina hiseq 2000          | 2,781,769  | Y          |
| nextseq 500                  | 2,320,288  | Y          |
| illumina hiseq 4000          | 1,701,041  | N          |
| nextseq 550                  | 1,359,672  | N          |
| hiseq x ten                  | 1,192,399  | Y          |
| nextseq 2000                 | 678,980    | N          |
| illumina hiseq 3000          | 383,170    | N          |
| illumina c                   | 205,644    | Y          |
| illumina novaseq x           | 172,756    | N          |
| nextseq 1000                 | 170,748    | N          |
| bgiseq-500                   | 152,831    | N          |
| illumina genome analyzer ii  | 152,149    | Y          |
| illumina novaseq x plus      | 146,379    | N          |
| illumina hiseq x             | 140,593    | Y          |
| dnbseq-g400                  | 124,940    | N          |
| dnbseq-t7                    | 119,742    | N          |
| illumina hiseq 1500          | 102,793    | N          |
| illumina genome analyzer iix | 96,622     | N          |
| mgiseq-2000rs                | 89,580     | N          |
| illumina iseq 100            | 62,647     | N          |
| illumina hiseq 1000          | 60,101     | Y          |
| illumina genome analyzer     | 53,061     | Y          |
| hiseq x five                 | 29,912     | Y          |
| dnbseq-g50                   | 6,544      | N          |
| bgiseq-50                    | 2,537      | N          |
| dnbseq-g400 fast             | 1,375      | N          |
| dnbseq-t10x4rs               | 600        | N          |
| onso                         | 237        | N          |

NOTE: Here we do not distinguish `hiseq x five`, `hiseq x ten` and `illumina hiseq x`.

## Collecting of read length data

The read length data are collected in batch.

Execute the following Python scripts in order:

```shell
python fetch_samples.py # Fetch sample metadata from ENA
python jmp_query.py   # Subsample and filtering
python fetch_rlen.py  # Fetch read length statistics
```

## Generating error profile files

### Illumina NovaSeq 6000

Up to 250bp PE.

```shell
mkdir -p NovaSeq6000_250bp
fasterq-dump --split-files --skip-technical --progress --threads 8 -A SRR30113744 -O NovaSeq6000_250bp/
../../opt/build_release_install/bin/art_profile_builder \
    --parallel 8 \
    --i-file NovaSeq6000_250bp/SRR30113744_1.fastq \
    --read_len 250 \
    --o-file1 NovaSeq6000_250bp/EmpNovaSeq6000_250R1.txt
../../opt/build_release_install/bin/art_profile_builder \
    --parallel 8 \
    --i-file NovaSeq6000_250bp/SRR30113744_2.fastq \
    --read_len 250 \
    --o-file1 NovaSeq6000_250bp/EmpNovaSeq6000_250R2.txt
```
