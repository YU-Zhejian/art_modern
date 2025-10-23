"""
This re-implementation of test-build.sh does not call make.
"""

from __future__ import annotations

import concurrent.futures
import os
import shutil
import subprocess
import tempfile
from typing import List, Optional
import threading
import argparse

MPIEXEC = os.environ.get("MPIEXEC")
MAKE = shutil.which("make")
CMAKE = shutil.which("cmake")
CTEST = shutil.which("ctest")
BASH = shutil.which("bash")
CMAKE_GENERATOR = "Ninja" if shutil.which("ninja") is not None else "Unix Makefiles"
SHDIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.environ.get("TEST_BUILD_LOG_DIR", os.path.abspath("test-build.log.d"))

IO_MUTEX = threading.Lock()


def probe_mkl() -> bool:
    # Dummy implementation for MKL probing
    with tempfile.TemporaryDirectory() as proj_dir, tempfile.TemporaryDirectory() as build_dir:
        shutil.copy(
            os.path.join(SHDIR, "test-build.d", "test-mkl.cmake"),
            os.path.join(proj_dir, "CMakeLists.txt"),
        )
        result = subprocess.run(
            [CMAKE, os.path.join(proj_dir, "CMakeLists.txt")],
            cwd=build_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return result.returncode == 0


def probe_concurrent_queue() -> Optional[str]:
    """
    Probes for the presence of the concurrentqueue library.

    :return: The path to the concurrentqueue header files if found, else None.
    """
    for p in [
        "/usr/include/concurrentqueue/moodycamel/concurrentqueue.h",
        "/usr/include/concurrentqueue/concurrentqueue.h",
    ]:
        if os.path.exists(p):
            return os.path.dirname(p)
    return None


class BuildConfig:
    def __init__(self, old_cmake_flags: List[str]):
        self.USE_RANDOM_GENERATOR = "PCG"
        self.USE_QUAL_GEN = "WALKER"
        self.USE_MALLOC = "AUTO"
        self.USE_THREAD_PARALLEL = "ASIO"
        self.USE_LIBFMT = "UNSET"
        self.USE_HTSLIB = "UNSET"
        self.USE_CONCURRENT_QUEUE = "UNSET"
        self.USE_ABSL = "UNSET"
        self.HELP_VERSION_ONLY = True
        self.old_cmake_flags = old_cmake_flags
        self.CMAKE_BUILD_TYPE = "Debug"

    def generate_cmake_opts(self) -> List[str]:
        cmake_flags = self.old_cmake_flags
        cmake_flags.append(f"-DCMAKE_BUILD_TYPE={self.CMAKE_BUILD_TYPE}")
        cmake_flags.append(f"-DUSE_RANDOM_GENERATOR={self.USE_RANDOM_GENERATOR}")
        cmake_flags.append(f"-DUSE_QUAL_GEN={self.USE_QUAL_GEN}")
        cmake_flags.append(f"-DUSE_MALLOC={self.USE_MALLOC}")
        cmake_flags.append(f"-DUSE_THREAD_PARALLEL={self.USE_THREAD_PARALLEL}")
        if self.USE_LIBFMT != "UNSET":
            cmake_flags.append(f"-DUSE_LIBFMT={self.USE_LIBFMT}")
        if self.USE_HTSLIB != "UNSET":
            cmake_flags.append(f"-DUSE_HTSLIB={self.USE_HTSLIB}")
        if self.USE_CONCURRENT_QUEUE != "UNSET":
            cmake_flags.append(f"-DUSE_CONCURRENT_QUEUE={self.USE_CONCURRENT_QUEUE}")
        if self.USE_ABSL != "UNSET":
            cmake_flags.append(f"-DUSE_ABSL={self.USE_ABSL}")
        return cmake_flags

    def copy(self) -> BuildConfig:
        new_config = BuildConfig(self.old_cmake_flags.copy())
        new_config.USE_RANDOM_GENERATOR = self.USE_RANDOM_GENERATOR
        new_config.USE_QUAL_GEN = self.USE_QUAL_GEN
        new_config.USE_MALLOC = self.USE_MALLOC
        new_config.USE_THREAD_PARALLEL = self.USE_THREAD_PARALLEL
        new_config.USE_LIBFMT = self.USE_LIBFMT
        new_config.USE_HTSLIB = self.USE_HTSLIB
        new_config.USE_CONCURRENT_QUEUE = self.USE_CONCURRENT_QUEUE
        new_config.USE_ABSL = self.USE_ABSL
        new_config.HELP_VERSION_ONLY = self.HELP_VERSION_ONLY
        return new_config


def do_build(config: BuildConfig, this_job_id: int) -> None:
    cmake_opts = config.generate_cmake_opts()
    build_dir = os.path.abspath(tempfile.mkdtemp(prefix="art_modern-test-build-"))
    install_dir = os.path.abspath(tempfile.mkdtemp(prefix="art_modern-test-build-install-"))
    with IO_MUTEX:
        print(f"{this_job_id} {' '.join(cmake_opts)}")
        print(f"{this_job_id} B: {build_dir}, I: {install_dir}")
    log_path = os.path.join(LOG_DIR, f"{this_job_id}.txt")
    with open(log_path, "wb") as log_file:
        log_file.write(f"Build log for job {this_job_id}\n".encode("utf-8"))
        log_file.flush()
        log_file.write(f"CMake options: {' '.join(cmake_opts)}\n".encode("utf-8"))
        log_file.flush()
        log_file.write(f"{this_job_id} B: {build_dir}, I: {install_dir}\n".encode("utf-8"))
        log_file.flush()
        # Invoke CMake to configure the build

        def run_wrapper(step_name: str, cmdline: List[str], *args, **kwargs):
            with IO_MUTEX:
                print(f"{this_job_id} {step_name} START")
            log_file.write(f"{this_job_id} {step_name} CMDLINE: { ' '.join(cmdline) }\n".encode("utf-8"))
            log_file.flush()
            log_file.write(f"{this_job_id} {step_name} START\n".encode("utf-8"))
            log_file.flush()
            if DRY_RUN:
                with IO_MUTEX:
                    print(f"{this_job_id} {step_name} DRYRUN")
                log_file.write(f"{this_job_id} {step_name} DRYRUN\n".encode("utf-8"))
                log_file.flush()
            else:
                proc = subprocess.run(cmdline, *args, **kwargs)
                if proc.returncode != 0:
                    with IO_MUTEX:
                        print(f"{this_job_id} {step_name} FAILED")
                    log_file.write(f"{this_job_id} {step_name} FAILED\n".encode("utf-8"))
                    log_file.flush()
                    raise RuntimeError(f"Step {step_name} failed for job {this_job_id}")
            with IO_MUTEX:
                print(f"{this_job_id} {step_name} DONE")
            log_file.write(f"{this_job_id} {step_name} DONE\n".encode("utf-8"))
            log_file.flush()

        run_wrapper(
            "CONFIG",
            [
                CMAKE,
                "-G",
                CMAKE_GENERATOR,
                "-Wdev",
                "-Wdeprecated",
                "--warn-uninitialized",
                "-DCEU_CM_SHOULD_ENABLE_TEST=ON",
                "-DCMAKE_VERBOSE_MAKEFILE=ON",
                *cmake_opts,
                f"-DCMAKE_INSTALL_PREFIX={install_dir}",
                os.getcwd(),
            ],
            cwd=build_dir,
            stdout=log_file,
            stderr=log_file,
        )

        run_wrapper(
            "BUILD",
            [
                CMAKE,
                "--build",
                build_dir,
                "-j",
                str(4),
            ],
            stdout=log_file,
            stderr=log_file,
        )
        run_wrapper(
            "CTEST",
            [CTEST, "--output-on-failure"],
            cwd=build_dir,
            stdout=log_file,
            stderr=log_file,
        )

        run_wrapper(
            "INSTALL",
            [CMAKE, "--install", build_dir],
            stdout=log_file,
            stderr=log_file,
        )

        run_wrapper(
            "TESTSMALL",
            [BASH, os.path.join(SHDIR, "test-small.sh")],
            stdout=log_file,
            stderr=log_file,
            env={
                **os.environ,
                "ART": os.path.join(install_dir, "bin", "art_modern" + ("-mpi" if WITH_MPI else "")),
                "HELP_VERSION_ONLY": "1" if config.HELP_VERSION_ONLY else "0",
                "MPIEXEC": MPIEXEC if WITH_MPI else "",
                "SAMTOOLS_THREADS": "4",
                "MPI_PARALLEL": "4",
                "PARALLEL": "2",
            },
        )

        shutil.rmtree(build_dir)
        log_file.write("RMDIR BUILD DONE\n".encode("utf-8"))
        log_file.flush()
        shutil.rmtree(install_dir)
        log_file.write("RMDIR INSTALL DONE\n".encode("utf-8"))
        log_file.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test build script re-implementation.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a dry run without executing commands.",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Use MPI for builds and tests.",
    )

    args, old_cmake_flags = parser.parse_known_args()
    DRY_RUN = args.dry_run
    WITH_MPI = args.mpi

    if WITH_MPI:
        if MPIEXEC is None:
            MPIEXEC = shutil.which("mpiexec")
        if MPIEXEC is None:
            raise RuntimeError("MPIEXEC not found in environment and mpiexec not in PATH")

    assert MAKE is not None, "make not found in PATH"
    assert CMAKE is not None, "cmake not found in PATH"
    assert CTEST is not None, "ctest not found in PATH"
    assert BASH is not None, "bash not found in PATH"

    if os.path.exists(LOG_DIR):
        shutil.rmtree(LOG_DIR)
    os.makedirs(LOG_DIR)
    RANDOM_GENERATORS = ["STL", "PCG", "BOOST"]

    print("Probing for MKL...")
    if probe_mkl():
        RANDOM_GENERATORS.append("MKL")

    print("Probing for concurrent queue...")
    concurrent_queue_path = probe_concurrent_queue()

    cmake_build_types = ["Debug", "Release", "RelWithDebInfo"]
    job_id = 0
    max_workers = min(5, os.cpu_count() // 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for cmake_build_type in cmake_build_types:
            bc = BuildConfig(old_cmake_flags)
            bc.CMAKE_BUILD_TYPE = cmake_build_type

            for use_thread_parallel in ["BS", "NOP"]:
                bc.USE_THREAD_PARALLEL = use_thread_parallel
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_THREAD_PARALLEL = "ASIO"

            for use_malloc in ["NOP", "JEMALLOC", "MIMALLOC"]:
                bc.USE_MALLOC = use_malloc
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_MALLOC = "AUTO"

            bc.USE_LIBFMT = "fmt"
            executor.submit(do_build, bc.copy(), job_id)
            job_id += 1
            bc.USE_LIBFMT = "UNSET"

            bc.USE_HTSLIB = "hts"
            executor.submit(do_build, bc.copy(), job_id)
            job_id += 1
            bc.USE_HTSLIB = "UNSET"

            if concurrent_queue_path is not None:
                bc.USE_CONCURRENT_QUEUE = concurrent_queue_path
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_CONCURRENT_QUEUE = "UNSET"

            bc.USE_ABSL = "SYS"
            executor.submit(do_build, bc.copy(), job_id)
            job_id += 1
            bc.USE_ABSL = "UNSET"

            bc.HELP_VERSION_ONLY = False

            for use_random_generator in RANDOM_GENERATORS:
                bc.USE_RANDOM_GENERATOR = use_random_generator
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_RANDOM_GENERATOR = "PCG"

            bc.USE_QUAL_GEN = "STL"
            executor.submit(do_build, bc.copy(), job_id)
            job_id += 1
            bc.USE_QUAL_GEN = "WALKER"
    executor.shutdown()
