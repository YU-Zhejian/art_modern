"""
Extract Gene Body Coverage Bias from 10x Genomics BAM files.

Algorithm:

1. Iterate over all reads in the BAM file
2. For each read, get the alignment position using TX tag
3. After processing all reads, calculate the effective gene body length for each transcript.
4. Calculate the gene body coverage bias for each transcript.

TODO: Get some file from the 10x Genomics website for viewing under the IGV.
"""
import pysam
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--out", required=True, help="Output file")
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    out = open(args.out, "w")
    for read in bam.fetch():
        if read.is_unmapped:
            continue
        # Required tags: TX (transcript), CB (cell barcode), UB (UMI barcode)
        if not read.has_tag("TX") or not read.has_tag("CB") or not read.has_tag("UB"):
            continue
        tx = read.get_tag("TX").split(";")[0]
        cb = read.get_tag("CB")
        ub = read.get_tag("UB")
        if not cb.endswith("-1"):
            continue
        cb = cb[:-2]
        # Get the alignment position
        transcript_id, match_pos, _ = tx.split(",")
        match_pos = int(match_pos)
        match_end = match_pos + read.query_alignment_length

