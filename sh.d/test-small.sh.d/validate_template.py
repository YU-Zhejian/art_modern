"""
Check whether the template BAM is valid.

Criteria:
- On SE: All reads start at contig start on pos strand.
- On PE: All first mates start at contig start on pos strand, all second mates start at contig end on neg strand.
- On MP: All first mates end at contig start on neg strand, all second mates end at contig end on pos strand.

TODO: Check C++ code and docs.
"""
