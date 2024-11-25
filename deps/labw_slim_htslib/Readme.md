# Readme for Modified `htslib`

Modifications done are:

- Removed GNU Autotools \& Makefile-based building system and use CMake instead.
- Dropped support for `libcurl` but enforce `libdeflate`, `libbz2` and `liblzma`. `libdeflate` acceleration details can be found at [official htslib benchmark](http://www.htslib.org/benchmarks/zlib.html).
- Removed plugin system.
- All files unused in `htscodecs` are removed.
- Support for external HTSCodecs and OpenSSL (For calculating MD5) are removed.

This modified `htsib` can process local BAM files at relative high speed. You may read more about `htslib` at [here](http://www.htslib.org/).
