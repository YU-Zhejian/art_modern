(am_compress-usage-section)=
# Usage of `am_compress`

`am_compress` is a compression tool that allows the user to compress arbitrary files to GZip and BGZip. Added in [1.3.4](#v-1.3.4-section).

## Implementation

This program in fact provides a user interface to `bgzf` data structure in HTSLib. Further compression algorithms can be added on request.

The input buffer size is always 4kb.

## Options

- `--o-compressed`: Output file name of the compressed data stream.
- `--o-compressed-compression`: Compression algorithm. Can be `gzip`, `bgzip` or `none`. If not set, will infer from extension.
- `--o-compressed-compression_level`: Compression level. Can be 1 (fastest) to 9 (best compression).
- `--o-compressed-buffer_size`: Buffer used for compression.
- `--o-compressed-num_threads`: Number of threads used in compression. Only used for BGZip algorithm.
- `--i-file`: Input file.

The program would emit warning if the output file extension does not match the desired algorithm.

## Usage Example

```shell
gzip --version
# gzip 1.12
# Copyright (C) 2018 Free Software Foundation, Inc.
# Copyright (C) 1993 Jean-loup Gailly.
# This is free software.  You may redistribute copies of it under the terms of
# the GNU General Public License <https://www.gnu.org/licenses/gpl.html>.
# There is NO WARRANTY, to the extent permitted by law.

# Written by Jean-loup Gailly.

libdeflate-gzip -V
# gzip compression program v1.24
# Copyright 2016 Eric Biggers

# This program is free software which may be modified and/or redistributed
# under the terms of the MIT license.  There is NO WARRANTY, to the extent
# permitted by law.  See the COPYING file for details.

pigz --version
# pigz 2.8

bgzip --version
# bgzip (htslib) 1.22.1
# Copyright (C) 2025 Genome Research Ltd.

igzip --version
# igzip command line interface 2.31.1

time gzip -k -f -9 SRR30113744_1.fastq
# real    5m21.201s
# user    5m20.821s
# sys     0m0.340s

time libdeflate-gzip -k -f -9 SRR30113744_1.fastq
# real    1m7.275s
# user    1m6.885s
# sys     0m0.292s

time pigz -k -f -9 SRR30113744_1.fastq
# real    0m16.392s
# user    5m19.205s
# sys     0m0.730s

time bgzip -k -f -l 9 -@ "$(nproc)" SRR30113744_1.fastq
# real    0m59.490s
# user    19m23.586s
# sys     0m1.468s

time igzip -k -f -3 -T "$(nproc)" SRR30113744_1.fastq
# real    0m6.299s
# user    0m4.569s
# sys     0m0.454s

time am_compress \
    --i-file SRR30113744_1.fastq \
    --o-compressed SRR30113744_1.fastq.gz \
    --o-compressed-compression_level 9 \
    --o-compressed-compression gzip \
    --o-compressed-num_threads 1
# With libdeflate as underlying implementation
# real    1m6.822s
# user    1m6.439s
# sys     0m0.347s

time am_compress \
    --i-file SRR30113744_1.fastq \
    --o-compressed SRR30113744_1.fastq.gz \
    --o-compressed-compression_level 9 \
    --o-compressed-compression gzip \
    --o-compressed-num_threads "$(nproc)"
# real    1m4.961s
# user    1m4.545s
# sys     0m0.412s

time am_compress \
    --i-file SRR30113744_1.fastq \
    --o-compressed SRR30113744_1.fastq.bgz \
    --o-compressed-compression_level 9 \
    --o-compressed-compression bgzip \
    --o-compressed-num_threads "$(nproc)"
# real    1m4.701s
# user    19m44.769s
# sys     0m22.885s
```
