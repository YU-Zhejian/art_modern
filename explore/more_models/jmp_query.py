import os

import numpy as np
import pandas as pd

MODELS = [
    "illumina novaseq 6000",
    "illumina hiseq 4000",
    "nextseq 550",
    "nextseq 2000",
    "illumina hiseq 3000",
    "illumina novaseq x",
    "nextseq 1000",
    "bgiseq-500",
    "illumina novaseq x plus",
    "dnbseq-g400",
    "dnbseq-t7",
    "illumina hiseq 1500",
    "illumina genome analyzer iix",
    "mgiseq-2000rs",
    "illumina iseq 100",
    "dnbseq-g50",
    "bgiseq-50",
    "dnbseq-g400 fast",
    "dnbseq-t10x4rs",
    "onso",
]


if __name__ == "__main__":
    dfs = []
    for model_name in MODELS:
        df = pd.read_csv(
            os.path.join("fetch_samples.d", f"{model_name.replace(' ', '_')}.tsv"),
            sep="\t",
            dtype={
                "run_accession": str,
                "fastq_bytes": str,
                "study_accession": str,
            },
        )
        if df.shape[0] == 0:
            print("  No records found, skipping.")
            continue
        num_records = df.shape[0]
        df = df.dropna()
        if df.shape[0] < 100:
            # Do not filter small datasets
            dfs.append(df)
            continue
        df["fastq_bytes_total"] = list(map(lambda x: sum([int(i) for i in x.split(";")]), df["fastq_bytes"]))
        df = (
            df.query(f"fastq_bytes_total >= ({1 * (1<<10)})")  # 1KB
            .query(f"fastq_bytes_total <= ({1 * (1<<30)})")  # 1GB
            .groupby("study_accession")
            .agg(np.random.choice)
            .reset_index(drop=True)
        )
        df["instrument_model"] = model_name
        print(
            f"{model_name}: Records after filtering: {num_records} -> {df.shape[0]} ({ df.shape[0] / num_records * 100 :.2f}%)"
        )
        dfs.append(df.sample(n=100) if df.shape[0] > 100 else df)
    df = pd.concat(dfs, ignore_index=True)
    df.to_parquet("filtered_samples.parquet", index=False)
    print(f"Total filtered samples: {df.shape[0]}")
