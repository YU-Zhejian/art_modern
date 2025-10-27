# Readme

The test would also check whether the produced FASTQ read ID is unique.

## Uniqueness Checking Algorithm

All FASTQ read IDs will firstly be separated based on their first 2 hex digits of CRC32 and compressed to GZip format.

## Dependencies

Requires zlib and libdeflate.

## TODO

A bloom filter is available [here](https://git.sr.ht/~sircmpwn/bloom/tree/master/item/main.c). However, we need to collect some hash functions.
