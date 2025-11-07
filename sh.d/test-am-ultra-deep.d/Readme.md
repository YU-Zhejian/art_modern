# Readme

The test would also check whether the produced FASTQ read ID is unique.

## Uniqueness Checking Algorithm

All FASTQ read IDs will firstly be separated based on their first 2 hex digits of CRC32 and compressed to GZip format.

**WARNING** Extensive memory usage.

## Dependencies

Requires zlib and libdeflate. Also requires `libxxhash-dev`.

## TODO

- Consider: <https://github.com/jwerle/murmurhash.c>
- Consider using FNV1a for key hashing algorithm in khash.
- Re-implement the bloom filter so those hash functions can be inlined.

## Used Third-Party Libraries

### Bloom Filter by Sam Hocevar

Available at [Git repo](https://git.sr.ht/~sircmpwn/bloom/) commit [`ee8657a4`](https://git.sr.ht/~sircmpwn/bloom/commit/ee8657a45522afaa5d6ca148c8ba1fa08220ccab).

Affected files:

- `bloom/*`

Under the Do What The Fuck You Want To Public License.

Some changes were made to the original code. The hash data type was changed to `uint32_t` and the hash function now requires length.

### RapidHash by Nicolas De Carli

Available at [Git repo](https://github.com/Nicoshev/rapidhash) commit [`2df1d03`](https://github.com/Nicoshev/rapidhash/commit/2df1d03aba25c4373f77dc0d559b2f789208fd42).

Affected files:

- `rapidhash/*`

Under the MIT license.

### KLib by Attractive Chaos

Available at [Git repo](https://github.com/attractivechaos/klib) commit [`5c1451c`](https://github.com/attractivechaos/klib/commit/5c1451caa1ee476624d00eed71810532c89b82d1).

Affected files:

- `klib/*`

Under the MIT license.

## Hedley by Evan Nemerson

Available at [Git repo](https://github.com/nemequ/hedley) commit [`8fb0604`](https://github.com/nemequ/hedley/commit/8fb0604a8095f6c907378cc3f0391520ae843f6f).

Affected files:

- `hedley/*`

Under CC0-1.0 license.
