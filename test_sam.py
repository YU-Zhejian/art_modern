"""
Test whether ART-generated SAM is correct.

Synopsis: $0 ref.fa file.bam
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
    BAM_CBACK = 9
    BAM_CBACK_STR = "B"

    STR_TO_INT = {
        BAM_CMATCH_STR: BAM_CMATCH,
        BAM_CINS_STR: BAM_CINS,
        BAM_CDEL_STR: BAM_CDEL,
        BAM_CREF_SKIP_STR: BAM_CREF_SKIP,
        BAM_CSOFT_CLIP_STR: BAM_CSOFT_CLIP,
        BAM_CHARD_CLIP_STR: BAM_CHARD_CLIP,
        BAM_CPAD_STR: BAM_CPAD,
        BAM_CEQUAL_STR: BAM_CEQUAL,
        BAM_CDIFF_STR: BAM_CDIFF,
        BAM_CBACK_STR: BAM_CBACK,
    }
    INT_TO_STR = {v: k for k, v in STR_TO_INT.items()}
    CONSUMES_QUERY = (True, True, False, False, True, False, False, True, True)
    CONSUMES_REFERENCE = (True, False, True, True, False, False, False, True, True)


if __name__ == "__main__":
    _, ref, alignment = sys.argv
    with pysam.FastaFile(ref) as ref_file, pysam.AlignmentFile(alignment, "rb") as alignment_file:
        alignments = list(alignment_file.fetch())
        if not alignments:
            raise ValueError("No alignments found")
        for aln in tqdm.tqdm(alignments):
            ref_seq = ref_file.fetch(aln.reference_name, aln.reference_start, aln.reference_end)
            query_seq = aln.query_sequence
            ref_ptr = 0
            query_ptr = 0
            for cigar_op, cigar_len in aln.cigartuples:
                if cigar_op == CigarOps.BAM_CEQUAL_STR:
                    assert ref_seq[ref_ptr : ref_ptr + cigar_len] == query_seq[query_ptr : query_ptr + cigar_len], aln
                elif cigar_op == CigarOps.BAM_CDIFF:
                    assert ref_seq[ref_ptr : ref_ptr + cigar_len] != query_seq[query_ptr : query_ptr + cigar_len], aln
                if CigarOps.CONSUMES_QUERY[cigar_op]:
                    query_ptr += cigar_len
                if CigarOps.CONSUMES_REFERENCE[cigar_op]:
                    ref_ptr += cigar_len
            assert ref_ptr == len(ref_seq), aln
            assert query_ptr == len(query_seq), aln
