# Readme for `art_modern`

Modernized ART that is parallelized and modularized using modern C++.

Changes:

- Supports 3 modes: `wgs`, `trans` and `templ`, similar to `pbsim3`.
- Supports 3 library construction methods: `se`, `pe` and `mp`.
- Masking detection dropped.
- Supports Illumina only.
- Support over the `aln` output format was dropped.
- Build system changed to CMake.
- All C++ code were re-implemented in C++14 with radical removal of duplicated or unused code.
- Random generator was changed from GNU Science Library (GSL) to Boost or C++ standard library.
- Logging re-implemented using Boost.
- Multithreading support implemented using Boost.
- Largely eliminated POSIX-only routines by Boost.
- Argument parser implemented in Boost.
- Built-in profiles are no longer supported. User must specify path to existing profiles.
- FASTA and SAM I/O re-implemented using `htslib` to allow huge FASTA files.

This simulator is based on the works of Weichun Huang <whduke@gmail.com>, under [GNU GPL v3](https://www.gnu.org/licenses/) license. The software is originally distributed [here](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) with following reference:

```bibtex
@article{10.1093/bioinformatics/btr708,
    author = {Huang, Weichun and Li, Leping and Myers, Jason R. and Marth, Gabor T.},
    title = "{ART: a next-generation sequencing read simulator}",
    journal = {Bioinformatics},
    volume = {28},
    number = {4},
    pages = {593-594},
    year = {2011},
    month = {12},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btr708},
    url = {https://doi.org/10.1093/bioinformatics/btr708},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/28/4/593/48879907/bioinformatics\_28\_4\_593.pdf},
}
```

Limitations:

TODO:

- Implement support over `sam` and `bam` output using `htslib`.
