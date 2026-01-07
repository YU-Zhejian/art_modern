"""
Usage: python validate_read_tags.py <input.sam> <expected tag 1> <expected tag 2> ...
"""

import sys
import pysam

if __name__ == "__main__":
    bam_name = sys.argv[1]
    if len(sys.argv) == 2:
        expected_tags = set()
    elif len(sys.argv) == 3:
        expected_tags = {sys.argv[2]}
    else:
        expected_tags = set(sys.argv[2:])

    with pysam.AlignmentFile(bam_name, "r", check_sq=False) as bam_file:
        for aln in bam_file:
            tags = set(dict(aln.get_tags()).keys())
            assert tags == expected_tags, f"{aln.query_name} has unexpected tags: {tags} != {expected_tags}"
