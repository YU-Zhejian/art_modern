import glob
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from get_profile import qual_max


def cmp(strand_no: int):
    profiles = glob.glob(f"profiles/*_{strand_no}.profile")
    profile_dict = {}
    for profile in profiles:
        profile_dict[profile.split("/")[-1].split(".")[0]] = pd.read_csv(profile, index_col=0)
    profile_1 = f"real_{strand_no}"
    for profile_2 in profile_dict.keys():
        for pos in range(0, 100):
            qual1 = profile_dict[profile_1].iloc[:, pos].to_numpy(dtype=float).ravel()
            qual2 = profile_dict[profile_2].iloc[:, pos].to_numpy(dtype=float).ravel()
            # mean1 = np.average(np.arange(0, qual_max), weights=qual1)
            mean2 = np.average(np.arange(0, qual_max), weights=qual2)
            correlation, p_value = mannwhitneyu(qual1, qual2)
            print(
                profile_1,
                profile_2.replace("_2", "").replace("_1", ""),
                pos,
                correlation,
                mean2,
                sep="\t",
            )


if __name__ == "__main__":
    print("profile_1", "profile_2", "pos", "correlation", "meanqual", sep="\t")
    cmp(1)
    cmp(2)
