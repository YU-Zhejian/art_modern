import glob
import json

import pandas as pd
import pyarrow as pa
import fastparquet as fp


if __name__ == "__main__":
    print("Versions: Pandas {}, PyArrow {}, Fastparquet {}".format(pd.__version__, pa.__version__, fp.__version__))

    all_files = glob.glob("monthly/*/*.parquet")
    odfs = []
    for fn in all_files:
        print(f"Processing {fn}...", end="")
        try:
            df = pd.read_parquet(fn, engine="pyarrow")
        except Exception as e:
            try:
                df = pd.read_parquet(fn, engine="fastparquet")
            except Exception as e:
                print("Failed")
                continue
        df_sel = df[["time", "data_source", "pkg_name", "counts"]]
        # Convert data_source to string to avoid issues with mixed types
        df_sel.loc[:, "data_source"] = df_sel.loc[:, "data_source"].astype(str)
        df_sel = df_sel.query("data_source == 'conda-forge' or data_source == 'bioconda'")
        # Sum the counts for each (time, data_source, pkg_name) triplet
        df_sel = df_sel.groupby(["time", "data_source", "pkg_name"], as_index=False, observed=False).sum("counts")
        odfs.append(df_sel)
        print("Done")
    if odfs:
        df_all = pd.concat(odfs, ignore_index=True)
        df_all.to_parquet("packing_status_all.parquet", index=False)
        print("Wrote packing_status_all.parquet with {} rows.".format(len(df_all)))
