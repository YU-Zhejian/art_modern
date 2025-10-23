"""
This re-implementation of test-build.sh does not call make.

TODO: Parallelization.
"""

import os
import shutil
import subprocess
import tempfile
from typing import List, Optional

OLD_CMAKE_FLAGS = os.environ.get("CMAKE_FLAGS")
MPIEXEC = os.environ.get("MPIEXEC")
MAKE = shutil.which("make")
CMAKE = shutil.which("cmake")
CTEST = shutil.which("ctest")
BASH = shutil.which("bash")
CMAKE_GENERATOR = "Ninja" if shutil.which("ninja") is not None else "Unix Makefiles"
WITH_MPI = True
SHDIR = os.path.dirname(os.path.abspath(__file__))


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
        self.BUILD_ONLY_TEST = 1
        self.old_cmake_flags = old_cmake_flags

    def generate_cmake_opts(self):
        cmake_flags = self.old_cmake_flags
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


def do_build(config: BuildConfig, log_path: str):
    cmake_opts = config.generate_cmake_opts()
    print("Building with the following CMake options: " + " ".join(cmake_opts))
    build_dir = os.path.abspath(tempfile.mkdtemp(prefix="art_modern-test-build-"))
    install_dir = os.path.abspath(
        tempfile.mkdtemp(prefix="art_modern-test-build-install-")
    )
    with open(log_path, "wb") as log_file:
        # Invoke CMake to configure the build
        subprocess.run(
            [
                CMAKE,
                "-G",
                CMAKE_GENERATOR,
                "-Wdev",
                "-Wdeprecated",
                "--warn-uninitialized",
                "-DCEU_CM_SHOULD_ENABLE_TEST=ON",
                *cmake_opts,
                f"-DCMAKE_INSTALL_PREFIX={install_dir}",
                os.getcwd(),
            ],
            cwd=build_dir,
            stdout=log_file,
            stderr=log_file,
            check=True,
        )
        # Build and test
        subprocess.run(
            [
                CMAKE,
                "--build",
                build_dir,
                "-j",
                20,  # TODO: make JOBS configurable
            ],
            stdout=log_file,
            stderr=log_file,
            check=True,
        )
        subprocess.run(
            [CTEST, "--output-on-failure"],
            cwd=build_dir,
            stdout=log_file,
            stderr=log_file,
            check=True,
        )
        subprocess.run(
            [CMAKE, "--install", build_dir],
            stdout=log_file,
            stderr=log_file,
            check=True,
        )
        subprocess.run(
            [BASH, os.path.join(SHDIR, "test-small.sh")],
            stdout=log_file,
            stderr=log_file,
            env = {
                **os.environ,
                "ART": os.path.join(install_dir, "bin", "art_modern" + ("-mpi" if MPIEXEC else "")),
                "HELP_VERSION_ONLY":  "0" if config.BUILD_ONLY_TEST == 1 else "1",
                "MPIEXEC": MPIEXEC if MPIEXEC else "",
            },
            check=True,
        )
    shutil.rmtree(build_dir)
    shutil.rmtree(install_dir)


if __name__ == "__main__":
    RANDOM_GENERATORS = ["STL", "PCG", "BOOST"]
    if probe_mkl():
        RANDOM_GENERATORS.append("MKL")
    cmake_build_types = ["Debug", "Release", "RelWithDebInfo"]
    for cmake_build_type in cmake_build_types:
        bc = BuildConfig([] ) # TODO
        bc.BUILD_ONLY_TEST = 1

        for use_thread_parallel in ["BS", "NOP"]:
            bc.USE_THREAD_PARALLEL = use_thread_parallel
            do_build(bc, f"build_log_{cmake_build_type}_{use_thread_parallel}.txt")
        bc.USE_THREAD_PARALLEL = "ASIO"

        for use_malloc in ["NOP", "JEMALLOC", "MIMALLOC"]:
            bc.USE_MALLOC = use_malloc
            do_build(bc, f"build_log_{cmake_build_type}_{use_malloc}.txt")
        bc.USE_MALLOC = "AUTO"

        bc.USE_LIBFMT = "fmt"
        do_build(bc, f"build_log_{cmake_build_type}_libfmt.txt")
        bc.USE_LIBFMT = "UNSET"

        bc.USE_HTSLIB = "hts"
        do_build(bc, f"build_log_{cmake_build_type}_htslib.txt")
        bc.USE_HTSLIB = "UNSET"

        bc.USE_CONCURRENT_QUEUE = probe_concurrent_queue()
        if bc.USE_CONCURRENT_QUEUE is not None:
            do_build(bc, f"build_log_{cmake_build_type}_concurrent_queue_path.txt")
        bc.USE_CONCURRENT_QUEUE = "UNSET"

        bc.USE_ABSL = "SYS"
        do_build(bc, f"build_log_{cmake_build_type}_absl.txt")
        bc.USE_ABSL = "UNSET"

        bc.BUILD_ONLY_TEST = 2

        for use_random_generator in RANDOM_GENERATORS:
            bc.USE_RANDOM_GENERATOR = use_random_generator
            do_build(bc, f"build_log_{cmake_build_type}_randgen_{use_random_generator}.txt")
        bc.USE_RANDOM_GENERATOR = "PCG"

        bc.USE_QUAL_GEN = "STL"
        do_build(bc, f"build_log_{cmake_build_type}_qualgen_stl.txt")
        bc.USE_QUAL_GEN = "WALKER"
