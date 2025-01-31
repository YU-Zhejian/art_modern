# Design Topics

Here lists various topics related to the design of `art_modern`.

## Parallelism

The parallelization strategy of different modes and input parsers are as follows:

| Parser \ Mode | `wgs`     | `trans`   | `templ`   |
|---------------|-----------|-----------|-----------|
| `memory`      | Coverage  | Batch     | Batch     |
| `htslib`      | Coverage  | **ERROR** | **ERROR** |
| `stream`      | **ERROR** | Batch     | Batch     |

Here, "coverage" means that all contigs are passed to all threads with coverage divided by the number of threads, while "batch" means that the entire data were split into batches (partitioning the existing in-memory data structure for `memory` parser or `--i-batch_size` for `stream` parser) and passed to different threads. The current batch-based design may be suboptimal if transcript/template lengths or coverages are not evenly distributed. You're recommended to shuffle the input data to avoid such problems using, i.e., `seqkit shuffle`. For even larger data, a proposed MPI-based parallelization strategy is in [TODO.md](TODO.md).

## Random Generators

### Bit Generation

The current random number generation function in each library is [MT19937](https://doi.org/10.1145/272991.272995), which may not be the best choice for performance-critical applications. However, it is the most widely used, well-known, be of moderate performance and cycle, and is implemented in all random number generator libraries (namely, Boost, GSL, STL, and Intel OneAPI MKL).

We choose not to use [Abseil Random](https://abseil.io/docs/cpp/guides/random) since (1) It is hard to bundle the entire Abseil random library with the project and (2) Its performance is not satisfying. Using Intel compiler, even absl::InsecureBitGen is slower than either boost::random::mt19937 or std::mt19937.

We choose not to use [cuRAND](https://docs.nvidia.com/cuda/curand/index.html) since it is hard to configure and its performance is not satisfying which may due to the time spent on copying data from GPU as the majority of our computation happens on CPU.

We may further introduce faster RNGs (e.g., PCG::pcg32) in the future. However, due to the fact that the current application is limited by IO, this is not a priority.

The current version does not support the specification of seed since it will result in reads with identical position and error profile in each thread. The current seed is set by the product of nanoseconds since epoch and hash of the current thread ID to allow each thread to generate different data.

### Distribution Sampling

The distribution sampling was currently implemented using each library's own implementation. For example; `boost::random::uniform_int_distribution` for Boost random generators, `std::uniform_int_distribution` for STL, `viRngUniform` for Intel MKL, and `gsl_rng_uniform_int` for GSL.

Further decoupling may be needed.

## I/O

The generated read data, which is represented as `PairwiseAlignment` class (`src/lib/PairwiseAlignment.hh`), will be passed to an output dispatcher which dispatches the data to all output writers registered to it.

Output writers are asynchronous. The `PairwiseAlignment`s are firstly formatted by each thread (since setting of `bam1_t*` is time-consuming) and then passed to a multi-producer single-consumer lock-free queue, which is then consumed by the file appender in a separate thread.

## Other Performance Bottlenecks

A majority of time was spent on generation of quality scores, which extensively calls `map::lower_bound()` function. Google B-tree map is used since it allows multiple values to be hold in one node, which improves performance when being compared red-black tree used in most STL implementations (including [EASTL](https://github.com/electronicarts/EASTL)).

## Building System

### Static Linking

The project should be able to be compiled into a fully static binary on [Alpine Linux](https://alpinelinux.org/) or [Void Linux](https://voidlinux.org/) with [musl libc](https://musl.libc.org/) as standard C library. See [this blog by Li Heng](https://lh3.github.io/2014/07/12/about-static-linking) for why static linking may simplify distribution and deployment of bioinformatics software.
