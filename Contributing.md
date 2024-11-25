# Readme for Developers

## Create Development Environment

Install [Conda](https://docs.conda.io/en/latest/) (Or [Mamba](https://mamba.readthedocs.io/en/latest/)/[MicroMamba](https://mamba.readthedocs.io/en/latest/micromamba.html)), and then execute:

```shell
conda env create -f art_modern.yml
```

## Parallelism

## Testing

Run `make testsmall` to run the tests using your default compiler.

## How to Submit Bugs \& Suggestions

### Performance Issues

If you are experiencing performance issues or if you're coming up with ideas for improving the performance of the code, please open an issue. You're recommended to report the following **in addition to** normal bug reports:

- The CMake commandline you used to build the code.
- A profiling report, such as from `valgrind`, `gprof`, `perf` or Intel Advisor.

