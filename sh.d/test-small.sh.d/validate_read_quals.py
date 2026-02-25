"""
Usage: python validate_read_tags.py <input.sam> <expected tag 1> <expected tag 2> ...
"""

import sys
import pysam

if __name__ == "__main__":
    bam_name = sys.argv[1]
    have_quals = sys.argv[2].lower() == "true"

    with pysam.AlignmentFile(bam_name, "r", check_sq=False) as bam_file:
        for aln in bam_file.fetch(until_eof=True):
            assert (aln.query_qualities is None) == (
                not have_quals
            ), f"{aln.query_name} has unexpected qualities presence: {aln.query_qualities is not None} != {have_quals}"
