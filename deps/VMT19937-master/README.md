# VMT19937

# A SIMD friendly pseudo random number generator, combining multiple Mersenne twister 19937 generators

## Introduction
This repository implements a new SIMD-friendly random number generator, which has the same statistical properties and period as the MT19937, the well known generator proposed by Nishimura and Matsumoto [^1].
It combines the random streams of multiple MT19937 instances with state vectors de-phased via jump-ahead transformations, then polls each instance in a round-robin fashion.
By evolving their vector states simultaneously, the new generator achieves perfect vectorization, fully leveraging on SIMD hardware capabilities.
Comprehensive test results demonstrate that the throughput of the new generator scales approximately linearly with the width of the SIMD registers used.

A paper describing in detail the implementation is available on [arXiv](https://arxiv.org/abs/2309.16682)

## Description
The VMT19937 generator is parametrized on 2 template parameters
```c++
template <size_t RegisterBitLen, VMT19937GenMode GenMode>
class VMT19937;
```

`RegisterBitLen` is the length of the registers used in bits.
It can be any of 32, 128, 256, 512. It reflects the capability of target architecture we target (e.g. if the target hardware only has support for SSE2 instructions, we can only use 32 and 128).
Performance improves with the register length. The set of available generators is determined automatically from compilation settings. For example, if we compile with the GCC compilation option _-mavx2_, then only 32, 128 and 256 are available. The compilation settings detection code is in the file [include/simd_config.h](include/simd_config.h).

`GenMode` determines the way in which the generator will be used. There are 3 modes to generate random numbers: one at a time, in blocks of 16 or in blocks with the same statesize as the generator. These modes correpond to the enum values below:
```c++
enum VMT19937GenMode { QM_Scalar, QM_Block16, QM_StateSize };
```
Performance improves with the length of the blocks. Usually `QM_Block16` is a good choice, as it represents a good compromize between performance and storage space.
The associated query functions are:
```c++
    // generates 1 random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32();

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // note the vector dst mus be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst);

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // note the vector dst mus be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst);
```

## Performance Stats
The table below shows the time in seconds to generate 5 billions of uniform discrete 32-bit random numbers in the range $[0,2^{32}-1]$.
`NBITS` and `GENMODE` are the template parameters of the generator. Peformance is compared against the original MT19937 generator and the vectorized SFMT19937 variation [^2].
```
-------------------------------------------------------------------
| GENERATOR | NBITS | GENMODE   | TARGET  | CPU-1 | CPU-2 | CPU-3 |
-------------------------------------------------------------------
| MT19937   | n.a.  | n.a.      | SSE2    | 31.56 | 20.07 | 16.90 | 
| SFMT19937 | n.a.  | n.a.      | SSE2    | 21.67 | 6.99  | 9.97  |
| VMT19937  | 32    | Scalar    | SSE2    | 20.83 | 11.10 | 13.54 |
| VMT19937  | 128   | Scalar    | SSE2    | 13.28 | 6.19  | 7.14  |
| VMT19937  | 128   | Block16   | SSE2    | 7.77  | 3.59  | 4.19  |
| VMT19937  | 128   | StateSize | SSE2    | 7.42  | 3.37  | 4.59  |
| VMT19937  | 256   | Scalar    | AVX     | n.a.  | 5.43  | 6.42  |
| VMT19937  | 256   | Block16   | AVX     | n.a.  | 2.15  | 2.15  |
| VMT19937  | 256   | StateSize | AVX     | n.a.  | 2.10  | 2.06  |
| VMT19937  | 512   | Scalar    | AVX512  | n.a.  | n.a.  | 5.66  |
| VMT19937  | 512   | Block16   | AVX512  | n.a.  | n.a.  | 1.45  |
| VMT19937  | 512   | StateSize | AVX512  | n.a.  | n.a.  | 1.14  |
-------------------------------------------------------------------
```
Performance stats obtained with the following CPUs:
- CPU-1: Intel(R) Celeron(R) J4125, cache 4Mb, frequency 2.0 GHz, burst frequency 2.7 GHz, SIMD support for SSE4.2 (a low end CPU).
- CPU-2: Intel® Core™ i9-12900H, cache 24Mb cache, base frequency 3.8GHz, turbo frequency 5.0 GHz, SIMD support for AVX2 (a high performance laptop CPU).
- CPU-3: Intel® Xeon® Gold 6234, cache 24.75Mb, base frequency 3.3GHz, turbo frequency 4.0 GHz, SIMD support for AVX512 (a high performance desktop CPU).

## Empirical random tests
The empirical quality of the pseudo random sequences generated via VMT19937 is tested with the _TestU01_ suite [^3]. The results are in the subfolder [logs/testu01](testu01-logs).
Results are comparable with the original MT19937 generator.

## Usage
The library is header only. You only need to include the header file [include/VMT19937.h](include/VMT19937.h).

Depending on your compilation settings, you have avaialble up to 12 generators, depending on the compination of the template parameters `RegisterBitLen` and `GenMode`.

You will need also some of the _jump-ahead_ matrices. These are available in _7z_ format in the [dat](dat) subfolder.
You can either uncompress the _7z_ files manually, or install _7za_ and then use the _make_ command.
For example, to extract the matrix _F19935.bits_, on Linux or Cygwin you can type the command `make dat/F19935.bits`.
Only one of the matrices is needed. Which one depends on the length of the SIMD registers available:
- 128-bits: _F19935.bits_
- 256-bits: _F19934.bits_
- 512-bits: _F19933.bits_

To initialize the generator, first we get the relevant matrix, then we initialize the generator.
For example, to instantiate a genrator with 128 bits regsiters and 
```c++
    // Get the relevant jump matrix. With 128 bits regsiter we pick F19935.bits
    const MT19937Matrix *jumpMatrix = new MT19937Matrix(std::string("./dat/F19935.bits"));

    // This is the initialization seed. Refer to the original MT19937 documentation.
    const uint32_t seedlength = 4;
    const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };
    
    // Create the generator setting QueryMode as Block16
    VMT19937<128, QM_Block16> mt(seedinit, seedlength, 0, nullptr, jumpMatrix);

    // The jump matrix is no longer needed and can be released here.
    delete jumpMatrix;
```

Then we can generater random numbers in blocks of 16
```c++
    // create storage vector aligned with cache lines, where we will store results
    AlignedVector<uint32_t, 64> buffer(16);

    // Query 10 times the generator in blocks of 16 numbers
    // We could also use mt.genrand_uint32(), which queries one number at a time, but it is slower.
    for (size_t i = 0; i < 10; ++i) {
        mt.genrand_uint32_blk16(buffer.data());
        for (size_t j = 0; j < 16; ++j)
            std::cout << buffer[j] << ", ";
    }
    std::cout << "\n";
```
The full source code for this example is in the file [src/demo.cpp](src/demo.cpp) in the routine `demo128`.

### Multiple independent generators
If you want to define multiple independent generator, for example to work with parallel Monte Carlo, you can use the common jump parameters of the constructor.
In that case you may need to extract further jump matrices, e.g. _F00100.bits_, which will generate independent streams with period of $2^{100}$.
```c++
    // Get the relevant jump matrix. With 128 bits registers we pick F19935.bits
    const MT19937Matrix *jumpMatrix = new MT19937Matrix(std::string("./dat/F19935.bits"));

    // Get the chosen common jump matrix. The size of the junp matrix determines the period
    // before the stream generated by multiple genarators become overlapping.
    // We can think about this as the period of the parallel generator.
    // In this example we use the F00100.bits matrix.
    const MT19937Matrix *commonJumpMatrix = new MT19937Matrix(std::string("./dat/F00100.bits"));

    // This is the initialization seed. Refer to the original MT19937 documentation.
    const uint32_t seedlength = 4;
    const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };
    
    // an array of parallel independent generators, which are guaraneteed not to be overlapping
    // up to a period of 2^100
    std::array<std::unique_ptr<VMT19937<128, QM_Block16>>, 10> parallelGenerators;

    // Create 10 multiple parallel generators with VecLen=128 and QueryMode=Block16
    for (size_t i = 0; i < 10; ++i)
        parallelGenerators[i].reset(new VMT19937<128, QM_Block16>(seedinit, seedlength, i, commonJumpMatrix, jumpMatrix));

    // The jump matrices are no longer needed and can be released here.
    delete jumpMatrix;
    delete commonJumpMatrix;
```

The full source code for this example is in the file [src/demo.cpp](src/demo.cpp) in the routine `demoParallel`.

## Requirements
The library is written in C++17. It uses extenseively Intel SIMD instructions, available on modern X86-64 processors.
It can be compiled for hardware with register length of 128, 256 or 512 bits.

## Build tests
Test and demo files are located in the [src](src) folder.
A [Makefile](Makefile) is provided for Linux or Cygwin. You need to define the `NBITS` environment variable. For example:
```
# compile with SSE2 code
make NBITS=128

# compile with AVX code
make NBITS=256

# compile with AVX512 code
make NBITS=512
```

## License
This is available with MIT [license](LICENSE). Note that the library includes also the [MT19937](mt19937-original) and [SFMT19937](SFMT-src-1.5.1) original source code.
These are distributed only for testing purpose and they have their own licenses.

## TestU01
To rerun the _TestU01_ tests, you need to download the _TestU01_ [package](http://simul.iro.umontreal.ca/testu01/tu01.html).
To build it, create an installation directory, e.g. _/workspace/repos/testu01/install_, then run the following commands:
```
./configure --prefix=/workspace/repos/testu01/install/ --disable-shared
make -j4
make install -j4
```
Then you need to build the VMT19937 test for _TESTU01_ and run the tests.
```
make -j4 NBITS=512 TESTU01_DIR=/workspace/repos/testu01/install/
PATH=/workspace/repos/testu01/install/bin:${PATH} make testu01logs
```

## References
[^1]: 1998, M. Matsumoto, T. Nishimura, _Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator_, ACM Trans. on Modeling and Computer Simulation, 8(1), 3-30. DOI: 10.1145/272991.272995.
[^2]: 2008, M. Saito and M. Matsumoto, _SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number Generator_, Monte Carlo and Quasi-Monte Carlo Methods 2006, Springer, 2008, pp. 607-622. DOI: 10.1007/978-3-540-74496-2\_36
[^3]: 2007, P. L'Ecuyer and R. Simard, _TestU01: A C Library for Empirical Testing of Random Number Generators_. ACM Transactions on Mathematical Software, Vol. 33, article 22.
