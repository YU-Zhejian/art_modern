import pyfastx

import sys

if __name__ == "__main__":
    in_fq_,targeted_rlen_str_ = sys.argv[1:]
    targeted_rlen_ = int(targeted_rlen_str_)
    for s in pyfastx.Fasta(in_fq_):
        if len(s) != targeted_rlen_:
            print(f"Contig {s.name} has length {len(s)}, expected {targeted_rlen_}")
            sys.exit(1)
    print("All contigs have the expected length.")
