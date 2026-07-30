# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.
- Add support for PE FASTA/FASTQ.

## Build and Test System


### `sh.d/test-build.py`

- `LOG_DIR` reads the environment variable `TEST_BUILD_LOG_DI` — a typo for `TEST_BUILD_LOG_DIR`. The override is therefore effectively unusable.
- `max_workers = min(5, os.cpu_count() // 4)` evaluates to `0` on machines with fewer than 4 CPUs, and `ThreadPoolExecutor(max_workers=0)` raises `ValueError`. Clamp to at least 1.
- The `run_ldd()` failure path writes nothing usable: `check_missing_libs()` logs to `PyLddSearch`'s own logger (`pyldd.log`) and raises a bare, message-less `AssertionError`, so the job's `ldd.log` is left 0 bytes and `ldd.lddtree.log` is never produced. The missing-library report belongs in the job log.
- `check_missing_libs()`: the `missing_libs = [...]` computation and the `raise` are nested one level too deep, inside the `for lib, libpath in libs_flattened` loop. The list is recomputed for every library and the function aborts on the first iteration instead of after logging all of them.
- Per-job `logging.FileHandler`s are never closed. On success the job directory is `rmtree`d while its handlers still hold the deleted files open, leaking file descriptors (3+ per job, tens of jobs per run).
- Timeouts (`BUILD` 360 s, `CTEST` 60 s, `INSTALL` 60 s, `TESTSMALL` 1200 s) and the build parallelism (`-j 4`) are hardcoded and unrelated to the machine. Combined with up to 5 concurrent workers this oversubscribes (5 x 4 = 20 compile jobs) *and* makes slower machines or LTO builds fail spuriously as `TIMEOUT`. These should be configurable and should honour the Makefile's `JOBS`.
- The `TESTSMALL` step runs `bash sh.d/test-small.sh` with a plain `os.environ`, unlike the Makefile's `testsmall-*` targets which go through `pixi run --frozen`. `samtools`, `htsfile`, `seqtk`, `fastqc`, `art-sam-validate` and the pixi `python` live only in `.pixi/envs/default/bin` and are not on `PATH` there, so `make testbuild-small` cannot work as written. Either re-exec under `pixi run`, or prepend the environment's `bin` to `PATH`.
- All jobs run `test-small.sh` concurrently against the *same* source tree: the same `PROJDIR/data/raw_data` (`gen_cov.py` regeneration) and the same exported `OPT_DIR`. This races between the up-to-5 workers.
- `BuildConfig.generate_cmake_opts()` does `cmake_flags = self.old_cmake_flags` — an alias, not a copy — and then appends to it, mutating the instance. Harmless today only because it is called exactly once per config; any retry would double-append every `-D` flag.
- In `run_ldd()`, `link_patterns = config.generate_link_patterns()` is assigned, discarded, and recomputed a few lines later; the `--dry-run` branch has a `FIXME` and consequently prints no patterns at all.
- `probe_using_cmake_script()` is dead code, and unlike `probe_using_cmake_project()` it never removes its log file on success.
- Cosmetic: `peobe_mkl_using_pkgconf` typo; `f"testmall.out.d"` should read `testsmall`; many `print(f"...")` have no placeholders; `LDD_CACHE` is a value assignment used as a type alias.
- `shutil.rmtree(LOG_DIR)` at startup deletes the directory unconditionally and without confirmation. Harmless with the PID-based default, dangerous once the `TEST_BUILD_LOG_DIR` override typo is fixed.
- `--fast-fail` raises from inside the worker; the handler in `__main__` calls `executor.shutdown(cancel_futures=True)` for *any* exception raised while iterating the futures inside the `with` block, and the reason surfaces only as `Exception inside TPE`.
- The generated matrix varies one option at a time against a fixed baseline for each of three build types, so the plain-default configuration is rebuilt many times while no cross-combination (e.g. `MIMALLOC` + `ONEMKL`) is ever tested.

### `Makefile`

- `clean` removes only `$(OPT_DIR) tmp build`. It leaves `cmake-build-clangd` (created by `lint`), the `test-build-*.log.d` directories that `test-build.py` leaves behind on failure, and `art_modern-log.d`. Since `pre-release-check.sh` begins with `make clean`, one stale log directory accumulates per failed run.
- `release-mpi` lacks `-DCEU_CM_SHOULD_ENABLE_TEST=ON` (and `-DBUILD_ART_MODERN_BENCHMARKS=ON`), both of which `release` has, so `testsmall-release-mpi` silently builds without tests. `debug` / `debug-mpi` do not have this asymmetry.
- `PIXI_PYTHON ?= $(shell pixi run --frozen which python)`: `?=` yields a recursively-expanded variable, so `pixi run` is re-executed on every single reference. Worse, if pixi is missing or the environment is not solved, the variable silently becomes empty and the recipe degrades to running `sh.d/test-build.py` as a bare command.
- `PACKAGE_VERSION ?= $(shell cat version.txt | tr -d '[:space:]')` has the same re-expansion problem, and it is `export`ed, so every sub-make re-runs it.
- `OPT_DIR ?= $(CURDIR)/opt` does not protect against `OPT_DIR=` being set-but-empty in the environment: `?=` will not override it, `mkdir -p $(OPT_DIR)/build_debug` then writes to `/build_debug`, and `clean` runs `rm -fr` over a bare ` tmp build`. Use `ifeq ($(strip $(OPT_DIR)),)`.
- `MAKEFLAGS += -R` does not do what its comment claims: built-in variables are defined before the makefile is parsed, so this only affects the sub-makes (`data/raw_data`, `docs/sphinx.d`, `packing`) — which may actually depend on them.
- `doc` and `serve-doc` invoke `pixi run --frozen -e doc python` while `help` uses `$(PIXI_PYTHON)` for the very same script. Pick one convention.
- The `testbuild*` targets pass `$(CMAKE_FLAGS)` as trailing positional arguments into `argparse.parse_known_args()`. Any flag beginning with `--` would be consumed or rejected by argparse instead of being forwarded to CMake.
- The `testbuild*` targets export `MPIEXEC` even for the non-MPI variants.
- There is no target that simply runs `ctest` against a `debug` / `release` build; unit tests are only ever exercised through `test-build.py`.

### `pre-release-check.sh`

- Every `make` invocation is suffixed with `|| true` and there is no summary or aggregate status, so **the script always exits 0** — a pre-release check that structurally cannot fail. Record per-stage results and exit non-zero.
- `TB_TARGET="testbuild" # -small`: the functional variant is commented out, so the pre-release check only ever compiles and runs `--help` / `--version`; the integration tests never run.
- No `set -u`; `. /opt/intel/oneapi/setvars.sh` is sourced with no existence check and no fallback, so on a machine without oneAPI the two Intel stages just fail into `|| true`.
- The two Intel stages are run through `run_gcc`, which is both a misnomer and the direct cause of the NCBI-NGS probe mismatch described above.
- `run_llvm` computes `CMAKE_TOOLCHAIN_FILE="$(pwd)/sh.d/toolchain/host-llvm/llvm-toolchain.cmake"`, so the script only works when invoked from the project root. It should `cd` to its own directory first.
- The hardcoded `~/opt/{sra-tools,ncbi-vdb}-3.3.0`, `~/opt/boost-1.89.0-clang` and `~/opt/fmt-12.0.0-clang` paths are never checked for existence, and the Boost version is duplicated between the `run_llvm` environment and the `-DBoost_DIR=` flag at the call sites.
- Blanket `# shellcheck disable=SC2317` at the top of the file hides genuinely unreachable code.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA, when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams, which allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worse if POSIX AIO is used.

## I/O Formats

- Support circular genome or RNA?
  - Search for new discussions in HTSLib, hts-specs, etc.
- Add `--i-nreads` to accurately specify the number of reads to simulate?

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
