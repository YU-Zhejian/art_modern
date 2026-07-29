import sys

import pysam

if __name__ == "__main__":
    in_fq_, targeted_rlen_str_ = sys.argv[1:]
    targeted_rlen_ = int(targeted_rlen_str_)
    with pysam.FastxFile(in_fq_) as fq_file:
        for entry in fq_file:
            seq_name = entry.name
            assert seq_name is not None, "Sequence name is None"
            seq = entry.sequence
            assert seq is not None, f"Sequence for {seq_name} is None"
            if len(seq) != targeted_rlen_:
                print(f"Contig {seq_name} has length {len(seq)}, expected {targeted_rlen_}")
                sys.exit(1)
    print("All contigs have the expected length.")
