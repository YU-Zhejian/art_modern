# Developing \& Contsibuting

## Common Development-Oriented Tasks

### Create Development Environment

Except [CMake](https://cmake.org) and other dependencies specified in [Install](Install.md), the development scripts also require other dependencies. Install [Conda](https://docs.conda.io/en/latest/) (Or [Mamba](https://mamba.readthedocs.io/en/latest/)/[MicroMamba](https://mamba.readthedocs.io/en/latest/micromamba.html)), and then execute:

```shell
conda env create -f art_modern.yml
```

### Testing

Some individual modules can be tested using CMake CTest system. Execute `ctest` in your building directory to test them.

Run `make testsmall` or `make testsmall-release` to run the integration tests using executables produced in `make debug` and `make release`.

### Profiling

We have `./profile.sh ${PROFILER}` to perform profiling of the software with profiler `${PROFILER}`. Profilers may include:

- `intel-advisor` for Intel Advisor. Intel compilers with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.
- `intel-vtune` for Intel VTune Profiler. Intel compilers with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.
- `nsys` for NVIDIA Nsight Systems. NVIDIA compiler with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.
- `valgrind` for Valgrind Callgrind. GCC with `Debug` mode will be used since optimization may produce instructions that are not supported by Valgrind.
- `amd-uprof` for AMD uProf. AMD compilers with `-mtune=native` and `-O3` with `RelWithDebInfo` mode will be used.

### Building Documentations

Create the development environment and make the documentation through:

```shell
conda activate art_modern
make doc
```

And the built documentations (HTML and PDF) should be in `doc/sphinx.d/_build/html` and `doc/sphinx.d/_build/latex` respectively. Note that for PDF output, you may need [latexmk](https://www.ctan.org/pkg/latexmk) and a working up-to-date [LaTeX](https://www.latex-project.org) distribution (e.g., [TeXLive](https://www.tug.org/texlive/), [MiKTeX](https://miktex.org/), [MacTeX](https://tug.org/mactex/)) installed.

### Others

- Run `make fmt` to format the code using [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html), [`sh`](https://github.com/mvdan/sh), [`cmake-format`](https://cmake-format.readthedocs.io/), and [`dos2unix`](https://www.freebsd.org/cgi/man.cgi?query=dos2unix&sektion=1).
- Run `make debug` to build the executable using debug mode with default compiler found by CMake.
- Run `make release` to build the executable using release mode with default compiler found by CMake.
- Run `make scc` to count lines of code written. Note that this excludes third-party codes so should be preferred over pure `scc` in project root.
- Run `make touch` to touch all files in the repository. This **MAY** work when CMake does strange things like compiling the source files again and again.
- Run `make raw_data` to download essential test data useful to various integration tests, benchmarks, etc.

## Get Engaged

### Issues

You're welcome to submit issues on GitHub if you've encountered any problems while using this software. You're recommended to take advantage of the bug report templates, which will help us to understand your problem better. You may also submit issues if you're unclear about the docs or have some ideas for improving the software itself.

### Pull Requests (PRs)

You're welcomed to send pull requests (PRs) to this project using the standard fork-and-pull-request workflow. However, please send an issue first to discuss the changes you're going to make (except spelling \& punctuation \& grammar issues). Otherwise, we may not be able to accept your PR.

Before you send a PR, please make sure that:

- `make testsmall` passes.
- You've run `make fmt` to format the code.
- You used Valgrind to check for memory leaks.

## Miscellaneous Developer-Oriented Documentation

```{toctree}
MakefileTargets.md
TODO.md
CODE_OF_CONDUCT.md
```
