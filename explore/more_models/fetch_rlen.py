import json
import os
import shutil
import subprocess

import pandas as pd
import tqdm

SRA_INFO_PARH = shutil.which("sra-info")
if SRA_INFO_PARH is None:
    raise ImportError("sra-info command not found in PATH. Please install SRA Toolkit.")


def get_sra_readlen(sra: str):
    # Use cmdline: sra-info -f json -S DRR298151
    cmdline = [SRA_INFO_PARH, "-f", "json", "-S", sra]
    s = subprocess.run(cmdline, capture_output=True, text=True, encoding="utf-8")
    if s.returncode != 0:
        raise RuntimeError(f"sra-info failed for {sra}: {s.stderr}")
    try:
        decoded = json.loads(s.stdout)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse JSON output for {sra}: {e}")
    read_lengths = decoded["SPOTS"]
    if len(read_lengths) != 1:
        raise ValueError(f"Expected one spot length, got {len(read_lengths)}")
    return [l["length"] for l in read_lengths[0]["reads"]]


if __name__ == "__main__":
    df = pd.read_parquet("filtered_samples.parquet")
    print(f"Total filtered samples: {df.shape[0]}")
    # Fetch read lengths
    if os.path.exists("filtered_samples_readlengths.parquet"):
        df_existing = pd.read_parquet("filtered_samples_readlengths.parquet")
        readlen_list = df_existing["read_lengths"].tolist()
    else:
        readlen_list = [[] for _ in range(df.shape[0])]
    for idx, row in tqdm.tqdm(df.iterrows(), total=df.shape[0]):
        sra = row["run_accession"]
        if readlen_list[idx] != []:
            continue  # Already fetched
        try:
            readlens = get_sra_readlen(sra)
            readlen_list[idx] = readlens
        except RuntimeError as e:
            print(f"  Failed to get read lengths for {sra}: {e}")
            # Do nothing, leave as empty list
        except ValueError as e:
            readlen_list[idx] = None
        if idx % 5 == 0:
            df_im = df.copy(True)
            df_im["read_lengths"] = readlen_list
            df_im.to_parquet("filtered_samples_readlengths.parquet", index=False)
    df["read_lengths"] = readlen_list
    df.to_parquet("filtered_samples_readlengths.parquet", index=False)
