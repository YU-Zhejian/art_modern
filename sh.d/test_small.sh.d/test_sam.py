"""
Test whether a SAM/BAM is correct.

Synopsis: $0 ref.fa file.bam

This script will work on sequences that are mapped to the opposite strand since in SAMv1.pdf, there is:

> For segments that have been mapped to the reverse strand, the recorded SEQ is reverse complemented from the original
> unmapped sequence and CIGAR, QUAL, and strand-sensitive optional fields are reversed and thus recorded
> consistently with the sequence bases as represented.

"""

import pysam
import sys
import tqdm


class CigarOps:
    BAM_CMATCH = 0
    BAM_CMATCH_STR = "M"
    BAM_CINS = 1
    BAM_CINS_STR = "I"
    BAM_CDEL = 2
    BAM_CDEL_STR = "D"
    BAM_CREF_SKIP = 3
    BAM_CREF_SKIP_STR = "N"
    BAM_CSOFT_CLIP = 4
    BAM_CSOFT_CLIP_STR = "S"
    BAM_CHARD_CLIP = 5
    BAM_CHARD_CLIP_STR = "H"
    BAM_CPAD = 6
    BAM_CPAD_STR = "P"
    BAM_CEQUAL = 7
    BAM_CEQUAL_STR = "="
    BAM_CDIFF = 8
    BAM_CDIFF_STR = "X"

    INT_TO_STR = [
        BAM_CMATCH_STR,
        BAM_CINS_STR,
        BAM_CDEL_STR,
        BAM_CREF_SKIP_STR,
        BAM_CSOFT_CLIP_STR,
        BAM_CHARD_CLIP_STR,
        BAM_CPAD_STR,
        BAM_CEQUAL_STR,
        BAM_CDIFF_STR,
    ]
    STR_TO_INT = {s: i for i, s in enumerate(INT_TO_STR)}
    CONSUMES_QUERY = (
        True,  # BAM_CMATCH_STR
        True,  # BAM_CINS_STR
        False,  # BAM_CDEL_STR
        False,  # BAM_CREF_SKIP_STR
        True,  # BAM_CSOFT_CLIP_STR
        False,  # BAM_CHARD_CLIP_STR
        False,  # BAM_CPAD_STR
        True,  # BAM_CEQUAL_STR
        True,  # BAM_CDIFF_STR
    )
    CONSUMES_REFERENCE = (
        True,  # BAM_CMATCH_STR
        False,  # BAM_CINS_STR
        True,  # BAM_CDEL_STR
        True,  # BAM_CREF_SKIP_STR
        False,  # BAM_CSOFT_CLIP_STR
        False,  # BAM_CHARD_CLIP_STR
        False,  # BAM_CPAD_STR
        True,  # BAM_CEQUAL_STR
        True,  # BAM_CDIFF_STR
    )


if __name__ == "__main__":
    _, ref, alignment = sys.argv
    flags = {
        "UNALIGNED": 0,
        "POS": 0,
        "NEG": 0,
    }
    with pysam.FastaFile(ref) as ref_file, pysam.AlignmentFile(alignment, "rb", check_sq=False) as alignment_file:
        alignments = list(alignment_file)
        if not alignments:
            raise ValueError("No alignments found")
        for aln in tqdm.tqdm(alignments):
            if aln.is_unmapped:
                flags["UNALIGNED"] += 1
                continue
            if aln.is_reverse:
                flags["NEG"] += 1
            else:
                flags["POS"] += 1
            ref_seq = ref_file.fetch(aln.reference_name, aln.reference_start, aln.reference_end).upper()
            query_seq = aln.query_sequence.upper()
            ref_ptr = 0
            query_ptr = 0
            genomic_ptr = aln.reference_start

            def where_we_are():
                return f"Q:{query_ptr}/R:{ref_ptr}/G:{genomic_ptr}/A:{aln.reference_name}:{aln.reference_start}- {aln.reference_end}:{'-' if aln.is_reverse else '+'}"

            for cigar_op, cigar_len in aln.cigartuples:
                if cigar_op == CigarOps.BAM_CEQUAL:
                    ref_seg = ref_seq[ref_ptr : ref_ptr + cigar_len]
                    query_seg = query_seq[query_ptr : query_ptr + cigar_len]
                    assert ref_seg == query_seg, f"{query_seg} != {ref_seg} {where_we_are()}"
                elif cigar_op == CigarOps.BAM_CDIFF:
                    ref_seg = ref_seq[ref_ptr : ref_ptr + cigar_len]
                    query_seg = query_seq[query_ptr : query_ptr + cigar_len]
                    assert ref_seg != query_seg, f"{query_seg} == {ref_seg} {where_we_are()}"
                if CigarOps.CONSUMES_QUERY[cigar_op]:
                    query_ptr += cigar_len
                if CigarOps.CONSUMES_REFERENCE[cigar_op]:
                    ref_ptr += cigar_len
                    genomic_ptr += cigar_len
            assert genomic_ptr == aln.reference_end, f"{genomic_ptr} != {aln.reference_end} {where_we_are()}"
            assert query_ptr == len(query_seq), f"{query_ptr} != {len(query_seq)} {where_we_are()}"
            assert ref_ptr == len(ref_seq), f"{ref_ptr} != {len(ref_seq)} {where_we_are()}"
    print(flags)
