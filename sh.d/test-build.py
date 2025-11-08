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
SUCCESS_ID = []
FAILED_ID = []


def probe_using_cmake_script(cmake_file_path: str) -> bool:
    # Dummy implementation for MKL probing
    with open(os.path.join(LOG_DIR, f"probe-{os.path.basename(cmake_file_path)}.log"), "wb") as log_file:
        result = subprocess.run(
            [CMAKE, "-P", os.path.join(SHDIR, "test-build.d", cmake_file_path)],
            stdout=log_file,
            stderr=log_file,
        )
    return result.returncode == 0


def probe_using_cmake_project(cmake_file_path: str) -> bool:
    # Dummy implementation for MKL probing
    with (
        tempfile.TemporaryDirectory() as proj_dir,
        tempfile.TemporaryDirectory() as build_dir,
        open(os.path.join(LOG_DIR, f"probe-{os.path.basename(cmake_file_path)}.log"), "wb") as log_file,
    ):
        shutil.copy(
            os.path.join(SHDIR, "test-build.d", cmake_file_path),
            os.path.join(proj_dir, "CMakeLists.txt"),
        )
        result = subprocess.run(
            [CMAKE, os.path.join(proj_dir, "CMakeLists.txt")],
            cwd=build_dir,
            stdout=log_file,
            stderr=log_file,
        )
        return result.returncode == 0


def probe_mkl() -> bool:
    return probe_using_cmake_project("test-mkl.cmake")


def probe_mimalloc() -> bool:
    return probe_using_cmake_project("test-mimalloc.cmake")


def probe_jemalloc() -> bool:
    return probe_using_cmake_script("test-jemalloc.cmake")


def probe_htshtslib() -> bool:
    return probe_using_cmake_script("test-htslib.cmake")


def probe_libfmt() -> bool:
    return probe_using_cmake_script("test-libfmt.cmake")


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

def probe_pcg_random_hpp() -> Optional[str]:
    """
    Probes for the presence of the pcg_random.hpp header file.

    :return: The path to the pcg_random.hpp header file if found, else None.
    """
    for p in [
        "/usr/include/pcg_random.hpp",
        "/usr/local/include/pcg_random.hpp",
    ]:
        if os.path.exists(p):
            return os.path.dirname(p)
    return None

class BuildConfig:
    def __init__(self, old_cmake_flags: List[str]):
        self.USE_RANDOM_GENERATOR = "PCG"
        self.USE_MALLOC = "AUTO"
        self.USE_THREAD_PARALLEL = "ASIO"
        self.USE_LIBFMT = "UNSET"
        self.USE_HTSLIB = "UNSET"
        self.USE_CONCURRENT_QUEUE = "UNSET"
        self.HELP_VERSION_ONLY = True
        self.old_cmake_flags = old_cmake_flags
        self.CMAKE_BUILD_TYPE = "Debug"
        self.AM_NO_Q_REVERSE = False

    def generate_cmake_opts(self) -> List[str]:
        cmake_flags = self.old_cmake_flags
        cmake_flags.append(f"-DCMAKE_BUILD_TYPE={self.CMAKE_BUILD_TYPE}")
        cmake_flags.append(f"-DUSE_RANDOM_GENERATOR={self.USE_RANDOM_GENERATOR}")
        cmake_flags.append(f"-DUSE_MALLOC={self.USE_MALLOC}")
        cmake_flags.append(f"-DUSE_THREAD_PARALLEL={self.USE_THREAD_PARALLEL}")
        if self.USE_LIBFMT != "UNSET":
            cmake_flags.append(f"-DUSE_LIBFMT={self.USE_LIBFMT}")
        if self.USE_HTSLIB != "UNSET":
            cmake_flags.append(f"-DUSE_HTSLIB={self.USE_HTSLIB}")
        if self.USE_CONCURRENT_QUEUE != "UNSET":
            cmake_flags.append(f"-DUSE_CONCURRENT_QUEUE={self.USE_CONCURRENT_QUEUE}")
        if self.AM_NO_Q_REVERSE:
            cmake_flags.append("-DAM_NO_Q_REVERSE=ON")
        if WITH_MPI:
            cmake_flags.append("-DWITH_MPI=ON")
        if self.AM_NO_Q_REVERSE:
            cmake_flags.append("-DAM_NO_Q_REVERSE=ON")
        return cmake_flags

    def copy(self) -> BuildConfig:
        new_config = BuildConfig(self.old_cmake_flags.copy())
        new_config.USE_RANDOM_GENERATOR = self.USE_RANDOM_GENERATOR
        new_config.USE_MALLOC = self.USE_MALLOC
        new_config.USE_THREAD_PARALLEL = self.USE_THREAD_PARALLEL
        new_config.USE_LIBFMT = self.USE_LIBFMT
        new_config.USE_HTSLIB = self.USE_HTSLIB
        new_config.USE_CONCURRENT_QUEUE = self.USE_CONCURRENT_QUEUE
        new_config.HELP_VERSION_ONLY = self.HELP_VERSION_ONLY
        new_config.AM_NO_Q_REVERSE = self.AM_NO_Q_REVERSE
        new_config.CMAKE_BUILD_TYPE = self.CMAKE_BUILD_TYPE
        return new_config


def do_build(config: BuildConfig, this_job_id: int) -> None:
    cmake_opts = config.generate_cmake_opts()
    build_dir = os.path.abspath(tempfile.mkdtemp(prefix="art_modern-test-build-"))
    install_dir = os.path.abspath(tempfile.mkdtemp(prefix="art_modern-test-build-install-"))
    with IO_MUTEX:
        print(f"{this_job_id} {' '.join(cmake_opts)}")
        print(f"{this_job_id} B: {build_dir}, I: {install_dir}")
    log_path = os.path.join(LOG_DIR, f"{this_job_id}.txt")
    with open(log_path, "wb", buffering=0) as log_file:
        log_file.write(f"Build log for job {this_job_id}\n".encode("utf-8"))
        log_file.write(f"CMake options: {' '.join(cmake_opts)}\n".encode("utf-8"))
        log_file.write(f"{this_job_id} B: {build_dir}, I: {install_dir}\n".encode("utf-8"))
        # Invoke CMake to configure the build

        def run_wrapper(step_name: str, cmdline: List[str], *args, **kwargs) -> bool:
            with IO_MUTEX:
                print(f"{this_job_id} {step_name} START")
            log_file.write(f"{this_job_id} {step_name} CMDLINE: { ' '.join(cmdline) }\n".encode("utf-8"))

            log_file.write(f"{this_job_id} {step_name} START\n".encode("utf-8"))

            if DRY_RUN:
                with IO_MUTEX:
                    print(f"{this_job_id} {step_name} DRYRUN")
                log_file.write(f"{this_job_id} {step_name} DRYRUN\n".encode("utf-8"))

            else:
                proc = subprocess.run(cmdline, *args, **kwargs)
                if proc.returncode != 0:
                    with IO_MUTEX:
                        print(f"{this_job_id} {step_name} FAILED")
                    log_file.write(f"{this_job_id} {step_name} FAILED\n".encode("utf-8"))
                    return False
            with IO_MUTEX:
                print(f"{this_job_id} {step_name} DONE")
            log_file.write(f"{this_job_id} {step_name} DONE\n".encode("utf-8"))
            return True

        if not run_wrapper(
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
        ):
            with IO_MUTEX:
                FAILED_ID.append(this_job_id)
            return

        if not run_wrapper(
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
        ):
            with IO_MUTEX:
                FAILED_ID.append(this_job_id)
                return
        if not run_wrapper(
            "CTEST",
            [CTEST, "--output-on-failure"],
            cwd=build_dir,
            stdout=log_file,
            stderr=log_file,
        ):
            with IO_MUTEX:
                FAILED_ID.append(this_job_id)
                return

        if not run_wrapper(
            "INSTALL",
            [CMAKE, "--install", build_dir],
            stdout=log_file,
            stderr=log_file,
        ):
            with IO_MUTEX:
                FAILED_ID.append(this_job_id)
                return

        if not run_wrapper(
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
        ):
            with IO_MUTEX:
                FAILED_ID.append(this_job_id)
                return

        shutil.rmtree(build_dir)
    os.remove(log_path)  # Remove log if successful.
    SUCCESS_ID.append(this_job_id)


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

    _args, _old_cmake_flags = parser.parse_known_args()
    DRY_RUN = _args.dry_run
    WITH_MPI = _args.mpi

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

    print("Probing for MKL...", end="")
    if probe_mkl():
        RANDOM_GENERATORS.append("ONEMKL")
        print("SUCCESS")
    else:
        print("FAIL")

    # TODO: Find SYSTEM_PCG.

    print("Probing for concurrent queue...", end="")
    concurrent_queue_path = probe_concurrent_queue()
    if concurrent_queue_path:
        print("SUCCESS")
    else:
        print("FAIL")

    MALLOCS = ["NOP"]
    print("Probing for jemalloc...", end="")
    if probe_jemalloc():
        MALLOCS.append("JEMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print("Probing for mimalloc...", end="")
    if probe_mimalloc():
        MALLOCS.append("MIMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for HTSLib...", end="")
    is_htslib_exist = probe_htshtslib()
    if is_htslib_exist:
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for {fmt}...", end="")
    is_libfmt_exist = probe_libfmt()
    if is_libfmt_exist:
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for <pcg_random.hpp>...", end="")
    pcg_random_hpp_exist_path = probe_pcg_random_hpp()
    if pcg_random_hpp_exist_path:
        print("SUCCESS")
    else:
        print("FAIL")

    cmake_build_types = ["Debug", "Release", "RelWithDebInfo"]
    job_id = 0
    max_workers = min(5, os.cpu_count() // 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for cmake_build_type in cmake_build_types:
            bc = BuildConfig(_old_cmake_flags)
            bc.CMAKE_BUILD_TYPE = cmake_build_type

            for use_thread_parallel in ["BS", "NOP"]:
                bc.USE_THREAD_PARALLEL = use_thread_parallel
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_THREAD_PARALLEL = "ASIO"

            for use_malloc in MALLOCS:
                bc.USE_MALLOC = use_malloc
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_MALLOC = "AUTO"

            if is_libfmt_exist:
                bc.USE_LIBFMT = "fmt"
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
                bc.USE_LIBFMT = "UNSET"

            if is_htslib_exist:
                bc.USE_HTSLIB = "hts"
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
                bc.USE_HTSLIB = "UNSET"

            if concurrent_queue_path is not None:
                bc.USE_CONCURRENT_QUEUE = concurrent_queue_path
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_CONCURRENT_QUEUE = "UNSET"

            bc.AM_NO_Q_REVERSE = True
            executor.submit(do_build, bc.copy(), job_id)
            job_id += 1
            bc.AM_NO_Q_REVERSE = False

            bc.HELP_VERSION_ONLY = False

            for use_random_generator in RANDOM_GENERATORS:
                bc.USE_RANDOM_GENERATOR = use_random_generator
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_RANDOM_GENERATOR = "PCG"
    executor.shutdown()
    print(f"Successful builds: {', '.join(SUCCESS_ID)}")
    print(f"Failed builds: {', '.join(FAILED_ID)}")
    if FAILED_ID:
        exit(1)
    else:
        exit(0)
