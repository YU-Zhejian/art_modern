import random
import os
import sys

import pyfastx

if __name__ == "__main__":
    rng = random.SystemRandom()
    basename = sys.argv[1]
    max_cov = int(sys.argv[2])
    with (
        open(f"{basename}.cov_strandless.tsv", "w", encoding="UTF-8") as strandless_cov_w,
        open(f"{basename}.cov_stranded.tsv", "w", encoding="UTF-8") as stranded_cov_w,
        open(f"{basename}.pbsim3.transcript", "w", encoding="UTF-8") as pbsim3_transcript_w,
    ):
        for name, seq in pyfastx.Fastx(f"{basename}.fa"):
            cov_pos = rng.random() * max_cov
            cov_neg = rng.random() * max_cov
            strandless_cov_w.write(f"{name}\t{cov_pos + cov_neg}\n")
            stranded_cov_w.write(f"{name}\t{cov_pos}\t{cov_neg}\n")
            pbsim3_transcript_w.write(f"{name}\t{cov_pos}\t{cov_neg}\t{seq}\n")
