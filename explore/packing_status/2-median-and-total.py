import pandas as pd
import pyarrow as pa
import fastparquet as fp


if __name__ == "__main__":
    print("Versions: Pandas {}, PyArrow {}, Fastparquet {}".format(
        pd.__version__, pa.__version__, fp.__version__
    ))
    df_all = pd.read_parquet("packing_status_all.parquet", engine="pyarrow")
    # Sum the median per (data_source, pkg_name) triplet
    df_median = df_all.groupby(
        ["data_source", "pkg_name"], as_index=False, observed=False
    )["counts"].median()

    df_median.to_parquet("packing_status_median.parquet", index=False)
    print("Wrote packing_status_median.parquet with {} rows.".format(len(df_median)))

    df_total = df_all.groupby(
        ["data_source", "pkg_name"], as_index=False, observed=False
    )["counts"].sum()
    df_total.to_parquet("packing_status_total.parquet", index=False)
    print("Wrote packing_status_total.parquet with {} rows.".format(len(df_total)))
