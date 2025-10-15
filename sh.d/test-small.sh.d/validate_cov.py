"""
FIXME: Stranded coverage not taken into account
"""
import sys
from enum import IntEnum
from typing import Union, Mapping

import pyfastx


class FileType(IntEnum):
    CONST_COV = 0
    COV_TSV = 1
    PBSIM3_TRANSCRIPT = 2


def parse_fq(in_fq: str, in_fa: str, file_type: FileType) -> Mapping[str, float]:
    if file_type == FileType.PBSIM3_TRANSCRIPT:
        contig_lengths = {}
        with open(in_fa) as fin:
            for l in fin:
                ls = l.strip().split("\t")
                contig_lengths[ls[0]] = len(ls[3])
    else:
        contig_lengths = {s.name: len(s) for s in pyfastx.Fasta(in_fa)}
    contig_bases = { contig_name: 0.0 for contig_name in contig_lengths }
    for seq_name, seq in pyfastx.Fastx(in_fq):
        contig_name = seq_name[1:].split(":")[0]
        contig_bases[contig_name ] += len(seq)
    for contig_name in contig_bases:
        contig_bases[contig_name] = contig_bases[contig_name] / contig_lengths[contig_name]
    return contig_bases

def parse_design(in_fa: str, in_cov: Union[float, str], file_type: FileType) -> Mapping[str, float]:
    if file_type == FileType.CONST_COV:
        contig_lengths = {s.name: len(s) for s in pyfastx.Fasta(in_fa)}
        return { contig_name: in_cov for contig_name in contig_lengths }
    elif file_type == FileType.COV_TSV:
        design_cov = {}
        with open(in_cov) as fin:
            for line in fin:
                ls = line.strip().split("\t")
                if len(ls) == 2:
                    contig_name, cov = ls
                    cov  = float(cov)
                elif len(ls) == 3:
                    contig_name, cov_pos, cov_neg = ls
                    cov = (float(cov_pos) + float(cov_neg))
                design_cov[contig_name] =cov
        return design_cov
    elif file_type == FileType.PBSIM3_TRANSCRIPT:
        design_cov = {}
        with open(in_fa) as fin:
            for l in fin:
                ls = l.strip().split("\t")
                design_cov[ls[0]] = float(ls[1]) + float(ls[2])
        return design_cov
    else:
        raise ValueError(f"Unknown file type: {file_type}")



def validate(in_fq: str, in_fa: str, in_cov: Union[float, str], file_type: FileType) -> None:
    fq_cov = parse_fq(in_fq, in_fa, file_type)
    design_cov = parse_design(in_fa, in_cov, file_type)
    for contig_name in fq_cov:
        if (fq_cov[contig_name] - design_cov[contig_name]) / design_cov[contig_name] > 0.02:
            raise ValueError(f"Coverage for {contig_name} is out of range: {fq_cov[contig_name]} vs {design_cov[contig_name]}")

if __name__ == "__main__":
    in_fq_, in_fa_, in_cov_, file_type_ = sys.argv[1:5]
    if file_type_ == "CONST_COV":
        file_type_ = FileType.CONST_COV
        in_cov_ = float(in_cov_)
    elif file_type_ == "COV_TSV":
        file_type_ = FileType.COV_TSV
    elif file_type_ ==  "PBSIM3_TRANSCRIPT":
        file_type_ = FileType.PBSIM3_TRANSCRIPT
    else:
        raise ValueError(f"Unknown file type: {file_type_}")

    validate(in_fq_, in_fa_, in_cov_, file_type_)
