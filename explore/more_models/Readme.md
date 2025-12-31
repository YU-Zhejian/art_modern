# Readme

Here records scripts and notes about additional models for `art_modern`.

## Existing Models

Here lists existing Illumina models supported by the original ART, which is still available in `art_modern`.

| Name | Name (`art_modern`)   | R1 `txt` file name    | R2 `txt` file name    | RLen1 | RLen2 | Model Name        |
|------|-----------------------|-----------------------|-----------------------|-------|-------|-------------------|
| GA1  | GA1Recalibrated_36bp  | `EmpR36R1`            | `EmpR36R2`            | 36    | 36    | GenomeAnalyzer I  |
| GA1  | GA1Recalibrated_44bp  | `EmpR44R1`            | `EmpR44R2`            | 44    | 44    | GenomeAnalyzer I  |
| GA2  | GA2Recalibrated_50bp  | `EmpR50R1`            | `EmpR50R2`            | 50    | 50    | GenomeAnalyzer II |
| GA2  | GA2Recalibrated_75bp  | `EmpR75R1`            | `EmpR75R2`            | 75    | 75    | GenomeAnalyzer II |
| HS10 | HiSeq1000_100bp       | `Emp100R1`            | `Emp100R2`            | 100   | 100   | HiSeq 1000        | 
| MSv1 | MiSeq_250bp           | `EmpMiSeq250R1`       | `EmpMiSeq250R2`       | 250   | 250   | MiSeq             |
| HS20 | HiSeq2000_100bp       | `HiSeq2kL100R1`       | `HiSeq2kL100R2`       | 100   | 100   | HiSeq 2000        |
| HS25 | HiSeq2500_125bp       | `HiSeq2500L125R1`     | `HiSeq2500L125R2`     | 126   | 126   | HiSeq 2500        |
| HS25 | HiSeq2500_150bp       | `HiSeq2500L150R1`     | `HiSeq2500L150R2`     | 150   | 150   | HiSeq 2500        |
| HSXn | HiSeqX_PCR_Free_150bp | `HiSeqXPCRfreeL150R1` | `HiSeqXPCRfreeL150R2` | 151   | 151   | HiSeqX PCR free   |
| HSXt | HiSeqX_TruSeq_150bp   | `HiSeqXtruSeqL150R1`  | `HiSeqXtruSeqL150R2`  | 151   | 151   | HiSeqX TruSeq     |
| MinS | MiniSeq_TruSeq_50bp   | `MiniSeqTruSeqL50`    | N/A                   | 51    | N/A   | MiniSeq TruSeq    |
| MSv3 | MiSeq_v3_250bp        | `MiSeqv3L250R1`       | `MiSeqv3L250R2`       | 251   | 251   | MiSeq v3          |
| NS50 | NextSeq500_v2_75bp    | `NextSeq500v2L75R1`   | `NextSeq500v2L75R2`   | 76    | 76    | NextSeq500 v2     |
| N/A  | GA1_36bp              | `Emp36R1`             | `Emp36R2`             | 36    | 36    | GenomeAnalyzer I  |

Following profiles are only available in `art_modern` but not in the original ART as a built-in model. However, the
model error files are provided in the original ART.

| Name (`art_modern`)     | R1 `txt` file name      | R2 `txt` file name      | RLen1 | RLen2 | Model Name        |
|-------------------------|-------------------------|-------------------------|-------|-------|-------------------|
| GA1_44bp                | `Emp44R1`               | `Emp44R2`               | 44    | 44    | GenomeAnalyzer I  |
| GA2_50bp                | `Emp50R1`               | `Emp50R2`               | 53    | 50    | GenomeAnalyzer II |
| GA2_75bp                | `Emp75R1`               | `Emp75R2`               | 75    | 75    | GenomeAnalyzer II |
| HiSeq2500Filtered_150bp | `HiSeq2500L150R1filter` | `HiSeq2500L150R2filter` | 150   | 150   | HiSeq 2500        |

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
| illumina hiscansq            | 12,770     | N/A        |
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
