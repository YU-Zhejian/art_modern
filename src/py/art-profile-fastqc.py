#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Based on the original FastQC source code:
<https://github.com/s-andrews/FastQC> commit 29c8f8b

From file: <uk/ac/babraham/FastQC/Graphs/BaseGroup.java>

Under the GPL v3 license

 * Copyright Copyright 2010-17 Simon Andrews
 *
 *    This file is part of FastQC.
 *
 *    FastQC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    FastQC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with FastQC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

import sys
from typing import Tuple, List

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib is required for this script", file=sys.stderr)
    sys.exit(1)

qual_max = 45

import numpy as np


def weighted_percentile(quals, qual_counts, percentile):
    quals = np.array(quals)
    qual_counts = np.array(qual_counts)
    sorter = np.argsort(quals)
    quals = quals[sorter]
    qual_counts = qual_counts[sorter]
    cum_counts = np.cumsum(qual_counts)
    total = cum_counts[-1]
    idx = np.searchsorted(cum_counts, percentile / 100 * total)
    return quals[min(idx, len(quals) - 1)]


def get_linear_interval(length: int) -> int:
    """
    Returns a sensible interval grouping for a given length.
    The first 9 positions are treated individually, then finds a grouping value
    which gives a total set of groups below 75.
    """
    base_values = [2, 5, 10]
    multiplier = 1

    while True:
        for base in base_values:
            interval = base * multiplier
            group_count = 9 + ((length - 9) // interval)
            if (length - 9) % interval != 0:
                group_count += 1
            if group_count < 75:
                return interval
        multiplier *= 10
        if multiplier == 10000000:
            raise Exception(f"Couldn't find a sensible interval grouping for length '{length}'")


def make_ungrouped_groups(max_length: int) -> List[Tuple[int, int]]:
    groups = []
    starting_base = 1
    interval = 1
    while starting_base <= max_length:
        end_base = starting_base + interval - 1
        if end_base > max_length:
            end_base = max_length
        groups.append((starting_base - 1, end_base))
        starting_base += interval
    return groups


def make_linear_base_groups(max_length: int) -> List[Tuple[int, int]]:
    if max_length <= 75:
        return make_ungrouped_groups(max_length)

    interval = get_linear_interval(max_length)
    starting_base = 1
    groups = []

    while starting_base <= max_length:
        end_base = starting_base + interval - 1

        if starting_base < 10:
            end_base = starting_base

        if starting_base == 10 and interval > 10:
            end_base = interval - 1

        if end_base > max_length:
            end_base = max_length

        groups.append((starting_base - 1, end_base))

        if starting_base < 10:
            starting_base += 1
        elif starting_base == 10 and interval > 10:
            starting_base = interval
        else:
            starting_base += interval
    return groups


def boxplot_metadata(quals, qual_counts):
    q1 = weighted_percentile(quals, qual_counts, 25)
    med = weighted_percentile(quals, qual_counts, 50)
    q3 = weighted_percentile(quals, qual_counts, 75)
    mean = np.average(quals, weights=qual_counts)

    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr

    # Whiskers: min/max within fences
    whislo = weighted_percentile(quals, qual_counts, 10)
    whishi = weighted_percentile(quals, qual_counts, 90)

    return {"whislo": whislo, "whishi": whishi, "q1": q1, "med": med, "q3": q3, "mean": mean, "fliers": []}


if __name__ == "__main__":
    # TODO: Add --sep_flag
    boxplot_data = []
    means = []
    quals_list = []
    qual_counts_list = []
    with sys.stdin as f:
        for l in f:
            if l.startswith("."):
                next_l = next(f)
                if not next_l.startswith("."):
                    break
                quals_list.append(list(map(int, l.strip().split("\t")[2:])))
                accumulated_counts = list(map(int, next_l.strip().split("\t")[2:]))
                not_accumulated_counts = [accumulated_counts[0]]
                for i in range(1, len(accumulated_counts)):
                    not_accumulated_counts.append(accumulated_counts[i] - accumulated_counts[i - 1])

                qual_counts_list.append(not_accumulated_counts)

    read_len = len(quals_list)
    lingrp = make_linear_base_groups(read_len)  # 0-based incl. excl.
    for lingrp_id, (lingrp_idx_start, lingrp_idx_end) in enumerate(lingrp):
        combined_quals = []
        combined_qual_counts = []
        for pos in range(lingrp_idx_start, lingrp_idx_end):
            if pos >= read_len:
                break
            combined_quals.extend(quals_list[pos])
            combined_qual_counts.extend(qual_counts_list[pos])
        combined_dedup_quals = []
        combined_dedup_qual_counts = []
        # Merge same quality scores
        appeared_quals = set(combined_quals)
        for possible_qual in range(qual_max):
            if possible_qual in appeared_quals:
                combined_dedup_quals.append(possible_qual)
                combined_dedup_qual_counts.append(0)
                for idx, (qual, count) in enumerate(zip(combined_quals, combined_qual_counts)):
                    if qual == possible_qual:
                        combined_dedup_qual_counts[-1] += count
        boxplot_data.append(boxplot_metadata(combined_dedup_quals, combined_dedup_qual_counts))
        means.append(boxplot_data[-1]["mean"])

    fig, ax = plt.subplots()
    ax.bxp(boxplot_data, positions=range(len(boxplot_data)), showmeans=False)
    x_names = []
    for lingrp_id, (lingrp_idx_start, lingrp_idx_end) in enumerate(lingrp):
        if lingrp_idx_start >= read_len:
            break
        if lingrp_idx_end > read_len:
            lingrp_idx_end = read_len
        if lingrp_idx_start + 1 == lingrp_idx_end:
            x_names.append(f"{lingrp_idx_start+1}")
        else:
            x_names.append(f"{lingrp_idx_start+1}-{lingrp_idx_end}")
    ax.set_xticks(range(len(boxplot_data)), labels=x_names, rotation=90)

    # Connect means
    ax.plot(range(len(boxplot_data)), means, color="red", marker="o", linestyle="-", label="Mean")

    ax.legend()
    ax.set_ylim(0, qual_max)
    plt.show()
