"""
This re-implementation of test-build.sh does not call make.

TODO: Add timeout. Also add time measurement.
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

THIS_PID = os.getpid()

MPIEXEC = os.environ.get("MPIEXEC")
MAKE = shutil.which("make")
CMAKE = shutil.which("cmake")
CTEST = shutil.which("ctest")
BASH = shutil.which("bash")
CMAKE_GENERATOR = "Ninja" if shutil.which("ninja") is not None else "Unix Makefiles"
SHDIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.environ.get("TEST_BUILD_LOG_DIR", os.path.abspath(f"test-build-{THIS_PID}.log.d"))
PKG_CONFIG_PATH = shutil.which("pkg-config")
if PKG_CONFIG_PATH is None:
    PKG_CONFIG_PATH = shutil.which("pkgconf")
if PKG_CONFIG_PATH is None:
    # Except
    raise RuntimeError("pkg-config or pkgconf not found in PATH")

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
    log_path = os.path.join(LOG_DIR, f"probe-{os.path.basename(cmake_file_path)}.log")
    with (
        tempfile.TemporaryDirectory() as proj_dir,
        tempfile.TemporaryDirectory() as build_dir,
        open(log_path, "wb") as log_file,
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
    if result.returncode == 0:
        os.unlink(log_path)
    return result.returncode == 0


def probe_pkg_config(package_name: str) -> bool:
    result = subprocess.run(
        [PKG_CONFIG_PATH, "--exists", package_name],
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


def peobe_mkl_using_pkgconf() -> List[str]:
    retl = []
    for p in ["mkl-sdl", "mkl-sdl-lp64"]:
        if probe_pkg_config(p):
            retl.append(p)
    return retl


class BuildConfig:
    def __init__(self, old_cmake_flags: List[str]):
        self.USE_RANDOM_GENERATOR = "PCG"
        self.USE_MALLOC = "AUTO"
        self.USE_THREAD_PARALLEL = "ASIO"
        self.USE_LIBFMT = None
        self.USE_HTSLIB = None
        self.USE_CONCURRENT_QUEUE = None
        self.HELP_VERSION_ONLY = True
        self.old_cmake_flags = old_cmake_flags
        self.CMAKE_BUILD_TYPE = "Debug"
        self.AM_NO_Q_REVERSE = False
        self.FIND_RANDOM_MKL_THROUGH_PKGCONF = None

    def generate_cmake_opts(self) -> List[str]:
        cmake_flags = self.old_cmake_flags
        cmake_flags.append(f"-DCMAKE_BUILD_TYPE={self.CMAKE_BUILD_TYPE}")
        cmake_flags.append(f"-DUSE_RANDOM_GENERATOR={self.USE_RANDOM_GENERATOR}")
        cmake_flags.append(f"-DUSE_MALLOC={self.USE_MALLOC}")
        cmake_flags.append(f"-DUSE_THREAD_PARALLEL={self.USE_THREAD_PARALLEL}")
        if self.USE_LIBFMT is not None:
            cmake_flags.append(f"-DUSE_LIBFMT={self.USE_LIBFMT}")
        if self.USE_HTSLIB is not None:
            cmake_flags.append(f"-DUSE_HTSLIB={self.USE_HTSLIB}")
        if self.USE_CONCURRENT_QUEUE is not None:
            cmake_flags.append(f"-DUSE_CONCURRENT_QUEUE={self.USE_CONCURRENT_QUEUE}")
        if self.AM_NO_Q_REVERSE:
            cmake_flags.append("-DAM_NO_Q_REVERSE=ON")
        if WITH_MPI:
            cmake_flags.append("-DWITH_MPI=ON")
        if self.AM_NO_Q_REVERSE:
            cmake_flags.append("-DAM_NO_Q_REVERSE=ON")
        if self.FIND_RANDOM_MKL_THROUGH_PKGCONF:
            cmake_flags.append(f"-DFIND_RANDOM_MKL_THROUGH_PKGCONF={self.FIND_RANDOM_MKL_THROUGH_PKGCONF}")
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
        new_config.FIND_RANDOM_MKL_THROUGH_PKGCONF = self.FIND_RANDOM_MKL_THROUGH_PKGCONF
        return new_config


def do_build(config: BuildConfig, this_job_id: int) -> None:
    cmake_opts = config.generate_cmake_opts()
    build_dir = os.path.join(LOG_DIR, f"art_modern-test-build-{this_job_id}")
    install_dir = os.path.join(LOG_DIR, f"art_modern-test-install-{this_job_id}")
    os.makedirs(build_dir)
    os.makedirs(install_dir)
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
                "ART_MODERN_PATH": os.path.join(install_dir, "bin", "art_modern" + ("-mpi" if WITH_MPI else "")),
                "APB_PATH": os.path.join(install_dir, "bin", "art_profile_builder" + ("-mpi" if WITH_MPI else "")),
                "HELP_VERSION_ONLY": "1" if config.HELP_VERSION_ONLY else "0",
                "MPIEXEC": MPIEXEC if WITH_MPI else "",
                "SAMTOOLS_THREADS": "4",
                "MPI_PARALLEL": "4",
                "PARALLEL": "2",
                "NO_FASTQC": "1",
                "OUT_DIR": os.path.join(LOG_DIR, f"test-small-out-{this_job_id}"),
            },
        ):
            with IO_MUTEX:
                FAILED_ID.append(this_job_id)
                return

        shutil.rmtree(build_dir)
        shutil.rmtree(install_dir)
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
    parser.add_argument(
        "--small",
        action="store_true",
        help="Run test_small.",
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
    print(f"LOG_DIR={LOG_DIR}")
    RANDOM_GENERATORS = ["STL", "PCG", "BOOST"]

    print("Probing for MKL through CMake...", end="")
    if probe_using_cmake_project("test-mkl.cmake"):
        RANDOM_GENERATORS.append("ONEMKL")
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for concurrent queue...", end="")
    concurrent_queue_path = probe_concurrent_queue()
    if concurrent_queue_path:
        print(concurrent_queue_path)
    else:
        print("FAIL")

    MALLOCS = ["NOP"]
    print("Probing for jemalloc...", end="")
    if probe_pkg_config("jemalloc"):
        MALLOCS.append("JEMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print("Probing for mimalloc...", end="")
    if probe_using_cmake_project("test-mimalloc.cmake"):
        MALLOCS.append("MIMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print("Probing for tcmalloc...", end="")
    if probe_pkg_config("libtcmalloc"):
        MALLOCS.append("TCMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print("Probing for tcmalloc_minimal...", end="")
    if probe_pkg_config("libtcmalloc_minimal"):
        MALLOCS.append("TCMALLOC_MINIMAL")
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for HTSLib...", end="")
    is_htslib_exist = probe_pkg_config("htslib")
    if is_htslib_exist:
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for {fmt}...", end="")
    is_libfmt_exist = probe_pkg_config("fmt")
    if is_libfmt_exist:
        print("SUCCESS")
    else:
        print("FAIL")

    print("Probing for <pcg_random.hpp>...", end="")
    pcg_random_hpp_exist_path = probe_pcg_random_hpp()
    if pcg_random_hpp_exist_path:
        RANDOM_GENERATORS.append("SYSTEM_PCG")
        print(pcg_random_hpp_exist_path)
    else:
        print("FAIL")

    print("Probing for MKL found through pkg-config...", end="")
    mkl_pc = peobe_mkl_using_pkgconf()
    if mkl_pc:
        print(";".join(mkl_pc))
    else:
        print("FAIL")

    cmake_build_types = ["Debug", "Release", "RelWithDebInfo"]
    job_id = 0
    max_workers = min(5, os.cpu_count() // 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for cmake_build_type in cmake_build_types:
            bc = BuildConfig(_old_cmake_flags)

            bc.HELP_VERSION_ONLY = not _args.small
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
                bc.USE_LIBFMT = None

            if is_htslib_exist:
                bc.USE_HTSLIB = "hts"
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
                bc.USE_HTSLIB = None

            if concurrent_queue_path is not None:
                bc.USE_CONCURRENT_QUEUE = concurrent_queue_path
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            bc.USE_CONCURRENT_QUEUE = None

            bc.AM_NO_Q_REVERSE = True
            executor.submit(do_build, bc.copy(), job_id)
            job_id += 1
            bc.AM_NO_Q_REVERSE = False

            for use_random_generator in RANDOM_GENERATORS:
                bc.USE_RANDOM_GENERATOR = use_random_generator
                executor.submit(do_build, bc.copy(), job_id)
                job_id += 1
            if mkl_pc:
                bc.USE_RANDOM_GENERATOR = "ONEMKL"
                for pc in mkl_pc:
                    bc.FIND_RANDOM_MKL_THROUGH_PKGCONF = pc
                    executor.submit(do_build, bc.copy(), job_id)
                    job_id += 1
                bc.FIND_RANDOM_MKL_THROUGH_PKGCONF = None

            bc.USE_RANDOM_GENERATOR = "PCG"
    executor.shutdown()
    SUCCESS_ID = sorted(SUCCESS_ID)
    FAILED_ID = sorted(FAILED_ID)
    print(f"Successful builds: {', '.join(map(str, SUCCESS_ID))}")
    print(f"Failed builds: {', '.join(map(str,FAILED_ID))}")
    if not FAILED_ID:
        shutil.rmtree(LOG_DIR)
    exit(1 if FAILED_ID else 0)
