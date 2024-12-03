# Design Topics

Here lists various topics related to the design of `art_modern`.

## Parallelism

The parallelization strategy of different modes and input parsers are as follows:

| Parser \ Mode | `wgs`     | `trans`   | `templ`   |
|---------------|-----------|-----------|-----------|
| `memory`      | Coverage  | Batch     | Batch     |
| `htslib`      | Coverage  | **ERROR** | **ERROR** |
| `stream`      | **ERROR** | Batch     | Batch     |

Here, "coverage" means that all contigs are passed to all threads with coverage divided by the number of threads, while "batch" means that the entire data were splitted into batches and passed to different threads. The current batch-based design may be sub-optimal if transcript/template lengths or coverages are not evenly distributed. You're recommended to shuffle the input data to avoid such problems using, i.e., `seqkit shuffle`. For even larger data, a proposed MPI-based parallelization strategy is in [TODO.md](TODO.md).

## Random Generators

The current random number generation function in each library is [MT19937](https://doi.org/10.1145/272991.272995), which may not be the best choice for performance-critical applications. However, it is the most widely used, well-known, be of moderate performance and cycle, and is implemented in all random number generator libraries (namely, Boost, GSL, STL, and Intel OneAPI MKL).

We may further introduce faster RNGs with shorter cycle (e.g., Taus) in the future for scientists that requires massive amount of data with less quality, and true rRNGs (e.g., `/dev/random`) for the contrary purpose. We may also further introduce cuRAND, which is said to be faster when generating a large amount of data. However, due to the fact that the current application is limited by IO, this is not a priority.

The current version does not support the specification of seed since it will result in reads with identical position and error profile in each thread. The current seed is set by the product of nanoseconds since epoch and hash of the current thread ID to allow each thread to generate different data.

## I/O

The generated read data, which is represented as `PairwiseAlignment` class (`src/lib/PairwiseAlignment.hh`), will be passed to an output dispatcher which dispatches the data to all output writers registered to it.

Output writers are asyncronous. The `PairwiseAlignment`s are firstly formatted by each thread (since setting of `bam1_t*` is time-consuming) and then passed to an multi-producer single-consumer lock-free queue, which is then consumed by the file appender in a separate thread.
