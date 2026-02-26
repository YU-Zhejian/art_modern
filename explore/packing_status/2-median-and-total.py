import datetime
import numpy as np
import pandas as pd
import pyarrow as pa
import fastparquet as fp


if __name__ == "__main__":
    print("Versions: Pandas {}, PyArrow {}, Fastparquet {}".format(pd.__version__, pa.__version__, fp.__version__))
    df_all = pd.read_parquet("packing_status_all.parquet", engine="pyarrow")

    df_median = df_all.groupby(["data_source", "pkg_name"], as_index=False, observed=False)["counts"].median()
    df_median.to_parquet("packing_status_median.parquet", index=False)
    print("Wrote packing_status_median.parquet with {} rows.".format(len(df_median)))

    df_total = df_all.groupby(["data_source", "pkg_name"], as_index=False, observed=False)["counts"].sum()
    df_total.to_parquet("packing_status_total.parquet", index=False)
    print("Wrote packing_status_total.parquet with {} rows.".format(len(df_total)))

    # Now calculating Conda-Velocity score.
    # Convert all months to year * 12 + month integer
    df_yr_mo_split = df_all["time"].str.split("-", expand=True)
    df_all["time_int"] = df_yr_mo_split[0].astype(int) * 12 + df_yr_mo_split[1].astype(int)
    now = datetime.datetime.now()
    current_time_int = now.year * 12 + now.month
    df_all["age_months"] = current_time_int - df_all["time_int"]
    # Filter all entries older than 12 months
    df_recent = df_all.query("age_months <= 12")
    # S = Sum_{mo} (counts_mo) * Exp[ -k * t]
    k = 0.15
    df_recent["decay_factor"] = (-k * df_recent["age_months"]).apply(lambda x: np.exp(x))
    df_recent["weighted_counts"] = df_recent["counts"] * df_recent["decay_factor"]
    df_velocity = df_recent.groupby(["data_source", "pkg_name"], as_index=False, observed=False)[
        "weighted_counts"
    ].sum()
    df_velocity = df_velocity.rename(columns={"weighted_counts": "counts"}).loc[
        :, ["data_source", "pkg_name", "counts"]
    ]
    df_velocity.to_parquet("packing_status_conda_velocity.parquet", index=False)
    print("Wrote packing_status_conda_velocity.parquet with {} rows.".format(len(df_velocity)))
