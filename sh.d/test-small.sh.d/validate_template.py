"""
Check whether the template BAM is valid.

Criteria:
- On SE: All reads start at contig start on pos strand.
- On PE: All first mates start at contig start on pos strand, all second mates start at contig end on neg strand.
- On MP: All first mates end at contig end on pos strand , all second mates end at contig start on neg strand.
"""

import sys
import pysam


def assess_se(aln: pysam.AlignedSegment, contig_length: int) -> bool:
    return not aln.is_reverse and aln.reference_start < 5 and aln.query_alignment_start < 5


def assess_pe(aln: pysam.AlignedSegment, contig_length: int) -> bool:
    if aln.is_read1:
        return not aln.is_reverse and aln.reference_start < 5 and aln.query_alignment_start < 5
    else:
        return (
            aln.is_reverse
            and abs(aln.reference_end - contig_length) < 5
            and abs(aln.query_alignment_end - len(aln.query_sequence)) < 5
        )


def assess_mp(aln: pysam.AlignedSegment, contig_length: int) -> bool:
    if aln.is_read1:
        return (
            not aln.is_reverse
            and abs(aln.reference_end - contig_length) < 5
            and abs(aln.query_alignment_end - len(aln.query_sequence)) < 5
        )
    else:
        return aln.is_reverse and aln.reference_start < 5 and aln.query_alignment_start < 5


if __name__ == "__main__":
    bam_name = sys.argv[1]
    mode = sys.argv[2].upper()  # SE, PE, MP
    with pysam.AlignmentFile(bam_name, "r", check_sq=False) as bam_file:
        for aln in bam_file:
            if aln.is_unmapped:
                print(f"Found unmapped read {aln.query_name}")
                sys.exit(1)
            contig_length = bam_file.get_reference_length(aln.reference_name)
            if mode == "SE":
                if not assess_se(aln, contig_length):
                    print(f"Read {aln.query_name} does not meet SE criteria")
                    sys.exit(1)
            elif mode == "PE":
                if not assess_pe(aln, contig_length):
                    print(f"Read {aln.query_name} does not meet PE criteria")
                    sys.exit(1)
            elif mode == "MP":
                if not assess_mp(aln, contig_length):
                    print(f"Read {aln.query_name} does not meet MP criteria")
                    sys.exit(1)
            else:
                print(f"Unknown mode: {mode}")
                sys.exit(1)
    print("All reads meet the criteria.")
