import glob
import json

import pandas as pd
import pyarrow as pa
import fastparquet as fp


if __name__ == "__main__":
    print("Versions: Pandas {}, PyArrow {}, Fastparquet {}".format(pd.__version__, pa.__version__, fp.__version__))

    all_files = glob.glob("monthly/*/*.parquet")
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
        schema = df.dtypes.to_dict()
        with open(fn.replace(".parquet", ".schema.json"), "w") as f:
            json.dump({k: str(v) for k, v in schema.items()}, f, indent=4)
            print(json.dumps({k: str(v) for k, v in schema.items()}))
        print(df.head(1))
