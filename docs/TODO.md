# TODO

## IMPORTANT

- For large-contigs test, the evenness of the coverage should be assessed.
- `make release` would fail on platforms without pkg-config, especially on Haiku OS and Debian GNU/Hurd.
- Add support for PE FASTA/FASTQ.

## Performance

- I/O:
  - The home-made "asynchronous IO" spent too much time in deallocating and creating new `std::unique_ptr`s. This problem is more obvious on smaller objects, like FASTA, when being compared to FASTQ.
  - Consider using the method implemented in `pigz`. That is, create a ring buffer that stores raw pointers to record datagrams, which allows reusing.
  - The current implementation passes too many small objects across the concurrent queue and I/O handlers, which is inefficient. This problem will be considerably worse if POSIX AIO is used.

## Restore the Documented Compiler Minimum by Dropping the Hard `<charconv>` Dependency

`docs/Install.md` claims GCC at least 7.4.0 and Clang at least 5.0.1, but the project cannot be built with either. `cigar_arr_to_str()` in `src/libam_support/utils/seq_utils.cc` calls `std::to_chars()` from `<charconv>`, which `libstdc++` only provides from 8.1; `libstdc++` 7 has no such header. This is core `libam_support`, not an optional component, so the build fails at the `#include`. Clang is affected identically, since on GNU/Linux it uses `libstdc++` by default and the standard library, not the compiler, sets the floor. Introduced together with the accelerated CIGAR-to-string conversion in [1.5.0](#v-1.5.0-section).

Either raise the documented minimum to GCC 9.1 (where C++17 library support is complete, and the same boundary as [the `<filesystem>` discussion](Design.md)), or keep the current minimum by making `<charconv>` optional:

- Guard the `#include` with `#if defined(__cpp_lib_to_chars) || __has_include(<charconv>)`, or set a macro from CMake, and provide a fallback that appends `std::to_string(len)` instead. That variant is already implemented and benchmarked as `cigar_arr_to_str_new` in `src/benchmark/cigar_to_str.cc`, so its cost relative to `std::to_chars()` is known rather than assumed.
- Keep both paths behind one function so that only `cigar_arr_to_str()` differs; nothing else in the project uses `<charconv>`.
- Whichever route is taken, add a compile check for the library features actually used. `deps/cmake_collections/test_cxx17.cmake` only probes generic C++17 and would not have caught this.

## Job Executor Lifetime Under the ASIO Backend

**Symptom:** with `USE_THREAD_PARALLEL=ASIO`, objects owned by a job appear not to be released once the job is done; `BS` does not show this.

**Source review (Boost 1.83, BS.thread\_pool 5.1.0) does not support that explanation.** Both backends destroy the task immediately after invoking it:

- ASIO. `post()` on `boost::asio::thread_pool` enqueues a `detail::executor_op` (`boost/asio/impl/thread_pool.hpp:137`). Its `do_complete()` moves the handler into a local, frees the operation's memory *before* the upcall, invokes, then destroys the local on return.
- BS. `worker()` does `task_t task = pop_task(); ... task();` inside an inner scope, so the `move_only_function` is destroyed at the closing brace.

So the `std::shared_ptr<JobExecutor>` captured by the lambda in `JobPool::add()` is released as soon as the task returns, under either backend.

**Real but bounded difference:** ASIO's `thread_info_base::deallocate()` stashes freed handler blocks in a per-thread `reusable_memory_[]` array instead of calling `free()`, releasing them only in `~thread_info_base()` (thread exit, i.e. `join()`). This is a fixed number of chunk-sized slots per thread, so it retains a constant, not something proportional to the number of jobs. It also only covers the operation shell; the handler destructor has already run.

**More likely cause — tasks that are destroyed without ever running.** `executor_op::do_complete()` skips the upcall when `owner == nullptr` (the shutdown path) and destroys the handler unrun. Such a job never enters `ArtJobExecutor::operator()()`, so `is_running_`, which is initialized to `true`, is still `true`. `JobPool::prune_finished_jobs()` keeps every executor for which `is_running()` holds, so those entries are never pruned, `cached_n_running_executors_` never reaches 0, and the debug check in `JobPool::stop()` fires. Externally this looks like leaked objects. Note also that `stop()` uses `pool_.wait()` for BS but `pool_.join()` for ASIO, which differ over queued-but-unstarted work.

**Next step:** replace `is_running_` with a tri-state (`NOT_STARTED` / `RUNNING` / `DONE`), or count entries into `operator()()`, and compare both backends. If ASIO shows posted-but-never-entered jobs, this is confirmed; if every job ran and the footprint still grows, the cause is elsewhere (look at the writer queues).

**Design note:** the lambda in `JobPool::add()` copies the `shared_ptr` while `executors_` holds another, so the pool is a genuine co-owner. Capturing a `std::weak_ptr`, or resetting the captured pointer at the end of the task, would make `prune_finished_jobs()` the sole release point and remove the backend dependence entirely.

## I/O Formats

- Support circular genome or RNA?
  - Search for new discussions in HTSLib, hts-specs, etc.
- Add `--i-nreads` to accurately specify the number of reads to simulate?

## Simulate Allele-Specific Expression

The simulator will accept a genome, a variation VCF file, and a gene annotation GTF file.
