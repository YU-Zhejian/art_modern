import sys

import numpy as np
import pandas as pd
import pyfastx
import tqdm

qual_max = 45
rlen = 100

if __name__ == "__main__":
    df = np.zeros((qual_max, rlen), dtype=int)
    fq = pyfastx.Fastx(sys.argv[1])
    for _, _, qual in tqdm.tqdm(fq):
        for pos in range(len(qual)):
            qual_int = max(min(ord(qual[pos]) - 33, qual_max - 1), 0)
            df[qual_int, pos] += 1
    pd.DataFrame(df, columns=range(rlen), index=range(qual_max)).to_csv(sys.argv[2])
