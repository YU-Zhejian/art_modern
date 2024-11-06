import random

import pyfastx

MRNA_PATH = "raw_data/ce11.mRNA_head.fa"
MAX_COV = 5

if __name__ == "__main__":
    rng = random.SystemRandom()
    with (
        open("raw_data/ce11.mRNA_head.cov_strandless.tsv", "w") as strandless_cov_w,
        open("raw_data/ce11.mRNA_head.cov_stranded.tsv", "w") as stranded_cov_w,
        open("raw_data/ce11.mRNA_head.pbsim3.transcript", "w") as pbsim3_transcript_w,
    ):
        for name, seq in pyfastx.Fastx(MRNA_PATH):
            cov_pos = rng.random() * MAX_COV
            cov_neg = rng.random() * MAX_COV
            strandless_cov_w.write(f"{name}\t{cov_pos + cov_neg}\n")
            stranded_cov_w.write(f"{name}\t{cov_pos}\t{cov_neg}\n")
            pbsim3_transcript_w.write(f"{name}\t{cov_pos}\t{cov_neg}\t{seq}\n")
