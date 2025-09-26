# README

David Blackman and Sebastiano Vigna introduced the [xoshiro/xoroshiro][] class of pseudorandom number generators, which are very efficient and, though they have relatively small states, display excellent statistical properties. The mathematical details are to be found in their [paper][].

`xoshiro.h` is a single header file, `C++` implementation of the _complete_ family of [xoshiro/xoroshiro][] pseudorandom uniform number generators.

The generators satisfy `C++`'s [`std::uniform_random_bit_generator`][] concept and can drive any distribution defined in the standard library.

We also provide efficient jump-ahead methods for *arbitrary* jump sizes, making these generators particularly suitable for parallel processing applications. They are excellent replacements for the "standard" _Mersenne Twister_ generator.

While the implementation is very general, there are simple type aliases for specific preferred instantiations that are known to work well in most situations.

## Example

Here is a simple example of using the _default xoshiro_ generator, `xso::rng`, to sample a normal distribution from the standard library.

We have copied the sample code from [`std::normal_distribution`][] webpage and used `xso::rng` as a drop-in replacement for [`std::mersenne_twister_engine`][].

```cpp
#include <xoshiro.h>
#include <random>
#include <map>
#include <iomanip>

int main()
{
    xso::rng gen;

    // lambda: draws a sample from a normal distribution and rounds it to an integer.
    std::normal_distribution d{5.0, 2.0};
    auto random_int = [&d, &gen] { return int(std::round(dist(gen))); };

    // Run many trials and create a histogram of the results in integer buckets.
    std::map<int, int> hist{};
    for (int n = 0; n != 25'000; ++n) ++hist[random_int()];

    // Print that histogram.
    for (auto [bucket, count] : hist) {
        auto n = static_cast<std::size_t>(count / 100);
        std::cout << std::setw(2) << bucket << ' ' << std::string(n, '*') << '\n';
    }
}
```

Everything is in the single `xoshiro.h` header file. The classes, functions etc., are all in the `xso` namespace. `xso::rng` is a type alias for `xso::rng64` which produces 64-bit outputs from a specific _xoshiro_ generator with 256 bits of state.

Here is the output from one run of the program:

```sh
-3
-2
-1
 0 **
 1 *******
 2 ****************
 3 *****************************
 4 ********************************************
 5 ************************************************
 6 ********************************************
 7 *****************************
 8 ****************
 9 *******
10 **
11
12
14
```

## Installation

This library is header-only, so there is nothing to compile & link—drop the `xoshiro.h` file somewhere convenient, and you are good to go.

Alternatively, if you are using `CMake`, you can use the standard `FetchContent` module by adding a few lines to your project's `CMakeLists.txt` file:

```cmake
include(FetchContent)
FetchContent_Declare(xoshiro URL https://github.com/nessan/bit/releases/download/current/xoshiro.zip)
FetchContent_MakeAvailable(xoshiro)
```

This command downloads and unpacks an archive of the current version of `xoshiro` to your project's build folder. You can then add a dependency on `xoshiro::xoshiro`, a `CMake` alias for `xoshiro`. `FetchContent` will automatically ensure the build system knows where to find the downloaded header files and any needed compiler flags.

Used like this, `FetchContent` will only download a minimal library version without any redundant test code, sample programs, documentation files, etc.

## Implementation

`C` versions of the generators are available on the author's [website][xoshiro/xoroshiro].
Wrapping those routines so they conform to the [`std::uniform_random_bit_generator`][] concept is a trivial exercise.

However, our implementation in `xoshiro.h` is distinguished in several other ways:

### Generality

Using `xoshiro.h`, you can create _any_ member of the _xoshiro/xoroshiro_ family.

We have `State` classes that are templatised across the number of state bits and the parameters (labelled `A`, `B`, and `C` in the literature) that determine how the state is advanced from step to step.

`Scrambler` classes are templatised across the other parameters (labelled `R`, `S`, and `T` in the literature) that determine how the higher dimensional state is scrambled/reduced to single 3-bit or 64-bit output words.

This means you can instantiate _any_ generator in the _xoshiro/xoroshiro_ family.

### No Compromise on Speed

For the reasonable optimisation levels you are likely to employ in any numerical code, the C++ versions perform identically to the simpler-looking `C` versions linked above.

### Simplicity

While you can instantiate _any_ generator in the _xoshiro/xoroshiro_ family, we also recognise that only a limited number of those generators have been vetted for suitability as being "good".

Therefore, we provide some _type aliases_ for the recommended default generators you should use in most cases.

### Arbitrary Jumps

We provide methods to advance a generator by _arbitrary_ and potentially vast numbers of steps. This contrasts with the `C` versions, which only define a limited number of jump possibilities.

Huge jumps are used to _partition_ a single random number stream into non-overlapping sub-streams. In parallel processing applications, the sub-streams drive independent jobs running on different compute cores.

### Sampling Methods

The `C++` standard library follows a typical design pattern for the facilities in its [`<random>`][] header. It maintains a strict separation between the classes that produce random bits and others that use those bits.

Uniform random bit generators produce streams of _uniformly_ distributed output words (usually 32-bit or 64-bit unsigned integers). Other classes and functions shape those core streams to simulate a desired distribution over some field of interest (for example, to generate variates from a uniform distribution of reals in the range 0 to 1).

The idea is reasonable enough. You can swap out the uniform random bit generator for a better one and continue to use the other functions without change.

However, this interface is quite complicated for the end user—particularly for a user who has no idea or interest in how all the generator/distribution machinery works.
Perhaps you google "`c++ random number generator`" and get advice to use the [`std::mersenne_twister_engine`][] class. However, that produces those unsigned output words that are useless by themselves. To simulate something as simple as a die roll, you must feed the thing into some other class and use that to get the required effect.

For this reason, we enrich our _xoshiro/xoroshiro_ classes with some utility methods that interface directly with the various distribution classes in the standard [`<random>`][] header.
That dice roll becomes something as simple as:

```cpp
xso::rng gen;
std::cout << "Six-sided dice roll: " << gen.roll() << '\n';
```

Similarly, you can ask one of our generators to flip a coin or shuffle the elements in a container.

Perhaps most importantly, the generators directly support the idea of _sampling_.
This includes pulling samples from a range, container, or an arbitrary distribution.

Here are some examples:

```cpp
std::cout << gen.sample(1, 10)   << '\n';   // <1>
std::cout << gen.sample(1., 10.) << '\n';   // <2>

std::normal_distribution nd{70., 15.};
std::cout << gen.sample(nd);                // <3>

std::array<double, 10> v;
gen.sample(nd, v.begin(), v.word_count());  // <4>

std::array<double, 5> u;
gen.sample(v, u.begin(), u.word_count());   // <5>

ge.shuffle(u);                              // <6>
```

1. Prints a random integer from $[1,10]$,  where each integer is equally likely to occur.
2. Prints a random real from a uniform distribution over $[1,10)$
3. Prints a random variate from a normal distribution with a mean of 70 and a standard deviation of 15.
4. Fills an array `v` with ten random variates from that same distribution.
5. Fills an array `u` with five elements drawn from `v` without replacement.
6. Shuffles the elements of `u`.

### Extra Analysis

Extra non-member functions for generator analysis are defined if the [`bit`][] library is available.

`bit` is a `C++` library for doing linear algebra over [GF(2)][] the simplest field of two elements $\{0,1\}$, where the usual arithmetic operations are performed mod 2.
The `bit` library is header only and is easily incorporated into any application.

If it is available, then `xoshiro.h` defines some extra functions that let you access the generator's _transition matrix_ and use/analyse it in various ways.

## Documentation

You can read the project's documentation [here](https://nessan.github.io/xoshiro/). \
The documentation site was generated using [Quarto](https://quarto.org).

### Contact

You can contact me by email [here](mailto:nzznfitz+gh@icloud.com).

### Copyright and License

Copyright (c) 2022-present Nessan Fitzmaurice. \
You can use this software under the [MIT license](https://opensource.org/license/mit).

<!-- Reference Links -->

[xoshiro/xoroshiro]: https://prng.di.unimi.it
[paper]: https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf
[`std::uniform_random_bit_generator`]: https://en.cppreference.com/w/cpp/numeric/random/uniform_random_bit_generator
[`std::normal_distribution`]: https://en.cppreference.com/w/cpp/numeric/random/normal_distribution
[`<random>`]: https://en.cppreference.com/w/cpp/header/random
[`std::mersenne_twister_engine`]: https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine
[`bit`]: https://nessan.github.io/bit
[GF(2)]: https://en.wikipedia.org/wiki/GF(2)
