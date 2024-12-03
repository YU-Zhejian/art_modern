# Readme for Modified `htslib`

Modifications done are:

- Removed GNU Autotools \& Makefile-based building system and use CMake instead.
- Dropped support for `libcurl`.
- Removed plugin system.
- All files unused in `htscodecs` are removed.
- Support for external HTSCodecs and OpenSSL (For calculating MD5) are removed.

This modified `htsib` can process local BAM files at relative high speed. You may read more about `htslib` at [here](http://www.htslib.org/).

## Bugs

The configuration system using CMake is still error-prone in Microsoft Windows.
