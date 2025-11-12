"""
FIXME: Stranded coverage not taken into account
"""

import statistics
import sys
from enum import IntEnum
from typing import Union, Mapping, Tuple

import pyfastx

import matplotlib.pyplot as plt


class FileType(IntEnum):
    CONST_COV = 0
    COV_TSV = 1
    PBSIM3_TRANSCRIPT = 2


def parse_fq(
    in_fq: str,
    in_fa: str,
    file_type: FileType,
    is_template: bool,
) -> Tuple[Mapping[str, float], Mapping[str, int], Mapping[str, int], Mapping[str, int]]:
    if file_type == FileType.PBSIM3_TRANSCRIPT:
        contig_lengths = {}
        with open(in_fa) as fin:
            for l in fin:
                ls = l.strip().split("\t")
                contig_lengths[ls[0]] = len(ls[3])
    else:
        contig_lengths = {s.name: len(s) for s in pyfastx.Fasta(in_fa)}
    if is_template:
        contig_lengths = {contig_name: 1 for contig_name in contig_lengths}
    contig_bases = {contig_name: 0 for contig_name in contig_lengths}
    contig_coverage = {contig_name: 0.0 for contig_name in contig_lengths}
    contig_nreads = {contig_name: 0 for contig_name in contig_lengths}
    for seq_name, seq, _ in pyfastx.Fastx(in_fq):
        contig_name = seq_name.split(":")[0]
        contig_bases[contig_name] += len(seq) if not is_template else 1
        contig_nreads[contig_name] += 1
    for contig_name in contig_bases:
        contig_coverage[contig_name] = contig_bases[contig_name] / contig_lengths[contig_name]
    return contig_coverage, contig_bases, contig_lengths, contig_nreads


def parse_design(in_fa: str, in_cov: Union[float, str], file_type: FileType) -> Mapping[str, float]:
    if file_type == FileType.CONST_COV:
        contig_lengths = {s.name: len(s) for s in pyfastx.Fasta(in_fa)}
        return {contig_name: in_cov for contig_name in contig_lengths}
    elif file_type == FileType.COV_TSV:
        design_cov = {}
        with open(in_cov) as fin:
            for line in fin:
                ls = line.strip().split("\t")
                if len(ls) == 2:
                    contig_name, cov = ls
                    cov = float(cov)
                elif len(ls) == 3:
                    contig_name, cov_pos, cov_neg = ls
                    cov = float(cov_pos) + float(cov_neg)
                design_cov[contig_name] = cov
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


def validate(in_fq: str, in_fa: str, in_cov: Union[float, str], file_type: FileType, is_template: bool) -> None:
    fq_cov, fq_bases, fq_contig_lengths, fq_contig_nreads = parse_fq(in_fq, in_fa, file_type, is_template)
    design_cov = parse_design(in_fa, in_cov, file_type)
    diffs = []
    diff_weights = []
    for contig_name in fq_cov:
        diff = (fq_cov[contig_name] - design_cov[contig_name]) / design_cov[contig_name]
        diffs.append(diff)
        diff_weights.append(design_cov[contig_name] * fq_contig_lengths[contig_name])
    # Geometric mean not used since the data is highly likely to contain zeros
    # stddev not used since the data may contain 1 data point only
    diff_mean = statistics.fmean(diffs, diff_weights)
    print("Diff: mean=", diff_mean, ".", sep="")
    if is_template:
        diff_threshold = 0.3
    else:
        diff_threshold = 0.1
    if (diff_mean) > diff_threshold:
        plt.scatter(fq_cov.values(), design_cov.values(), alpha=0.01)
        plt.show()
        plt.clf()
        raise ValueError("Mean coverage difference is greater than 10%.")


if __name__ == "__main__":
    in_fq_, in_fa_, in_cov_, file_type_, is_template_ = sys.argv[1:6]
    if file_type_ == "CONST_COV":
        file_type_ = FileType.CONST_COV
        in_cov_ = float(in_cov_)
    elif file_type_ == "COV_TSV":
        file_type_ = FileType.COV_TSV
    elif file_type_ == "PBSIM3_TRANSCRIPT":
        file_type_ = FileType.PBSIM3_TRANSCRIPT
    else:
        raise ValueError(f"Unknown file type: {file_type_}")

    validate(in_fq_, in_fa_, in_cov_, file_type_, is_template_ == "IS_TEMPLATE")
