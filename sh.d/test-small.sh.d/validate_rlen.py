import pyfastx

import sys

if __name__ == "__main__":
    in_fq_, targeted_rlen_str_ = sys.argv[1:]
    targeted_rlen_ = int(targeted_rlen_str_)
    for seq_name, seq, _ in pyfastx.Fastx(in_fq_):
        if len(seq) != targeted_rlen_:
            print(f"Contig {seq_name} has length {len(seq)}, expected {targeted_rlen_}")
            sys.exit(1)
    print("All contigs have the expected length.")
