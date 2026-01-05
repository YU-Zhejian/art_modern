"""
This re-implementation of test-build.sh does not call make.
"""

from __future__ import annotations

import concurrent.futures
import glob
import itertools
import os
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import argparse
import logging

from elftools.elf.elffile import ELFFile
from elftools.elf.dynamic import DynamicSection

from typing import List, Optional, Tuple, Dict

THIS_PID = os.getpid()

MPIEXEC = os.environ.get("MPIEXEC")
MAKE = shutil.which("make")
CMAKE = shutil.which("cmake")
CTEST = shutil.which("ctest")
BASH = shutil.which("bash")
CMAKE_GENERATOR = "Ninja" if shutil.which("ninja") is not None else "Unix Makefiles"
SHDIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.environ.get("TEST_BUILD_LOG_DI", os.path.abspath(f"test-build-{THIS_PID}.log.d"))
PKG_CONFIG_PATH = shutil.which("pkg-config")
if PKG_CONFIG_PATH is None:
    PKG_CONFIG_PATH = shutil.which("pkgconf")
if PKG_CONFIG_PATH is None:
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


def probe_ncbi_ngs() -> bool:
    p = subprocess.run(
        [BASH, os.path.join(SHDIR, "test-build.d", "have_ncbi_ngs.sh")],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return p.returncode == 0


class PyLddSearch:
    ld_so_paths: List[str]
    ld_library_path: List[str]
    ld_preload: List[str]

    lh: logging.Logger

    def _filter_ld_path(self, ld_paths: List[str]) -> List[str]:
        ld_paths_filtered = []
        for p in ld_paths:
            p = p.strip()
            if not p:  # Ignore blank ones
                continue
            if os.path.exists(p):
                abspath = os.path.abspath(p)
                if abspath not in ld_paths_filtered:
                    ld_paths_filtered.append(abspath)
                else:
                    self.lh.warning("LD path DUPLICATED: %s", abspath)
            else:
                self.lh.warning("LD path ENOENT: %s", p)
        return ld_paths_filtered

    def __init__(self, ld_so_conf_path: str = "/etc/ld.so.conf"):
        self.lh = logging.getLogger(self.__class__.__name__)
        self.lh.handlers = []
        self.lh.setLevel(logging.INFO)
        self.lh.addHandler(logging.StreamHandler(sys.stderr))
        self.lh.addHandler(logging.FileHandler(os.path.join(LOG_DIR, "pyldd.log"), "w", "utf-8"))

        for handler_ in self.lh.handlers:
            handler_.setLevel(logging.INFO)
            handler_.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
        # Read /etc/ld.so.conf and parse it
        self.ld_so_paths = self._filter_ld_path(self._parse_ld_so_conf(ld_so_conf_path))
        # Read LD_LIBRARY_PATH
        self.ld_library_path = self._filter_ld_path(os.environ.get("LD_LIBRARY_PATH", "").split(":"))
        self.lh.info("Final Resolved and valid /etc/ld.so.conf: %s", ":".join(self.ld_so_paths))
        self.lh.info("Final Resolved and valid LD_LIBRARY_PATH: %s", ":".join(self.ld_library_path))

        # Now parsing LD_PRELOAD
        ld_preload = os.environ.get("LD_PRELOAD", "")
        ld_preload_filtered = []
        if ld_preload:
            for p in ld_preload.split(":"):
                p = p.strip()
                if not p:
                    continue
                if os.path.exists(p):
                    abspath = os.path.abspath(p)
                    if abspath not in ld_preload_filtered:
                        ld_preload_filtered.append(abspath)
                    else:
                        self.lh.warning("LD preload DUPLICATED: %s", abspath)
                else:
                    # Try to search in ld_paths
                    found_path = self._search_library_in_paths(p, [], [])
                    if found_path is not None:
                        if found_path not in ld_preload_filtered:
                            ld_preload_filtered.append(found_path)
                        else:
                            self.lh.warning("LD preload DUPLICATED: %s", found_path)
                    else:
                        self.lh.warning("LD preload ENOENT: %s", p)
        self.ld_preload = ld_preload_filtered
        self.lh.info("Final Resolved and valid LD preloads (LD_PRELOAD): %s", ":".join(self.ld_preload))

        # We do NOT parse ld.so.cache.

    def _parse_ld_so_conf(self, ld_so_conf_path: str) -> List[str]:
        paths = []
        if not os.path.exists(ld_so_conf_path):
            return paths

        with open(ld_so_conf_path, "r", encoding="UTF-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if line.startswith("include "):
                    # Handle 'include /etc/ld.so.conf.d/*.conf'
                    pattern = line.split(maxsplit=1)[1]
                    for included_file in glob.glob(pattern):
                        paths.extend(self._parse_ld_so_conf(included_file))
                else:
                    paths.append(line)
        return paths

    def _search_library_in_paths(
        self,
        lib_name: str,
        dt_rpath: List[str],
        dt_runpath: List[str],
    ) -> Optional[str]:
        """
        Search for the given library name in the resolved ld_paths.
        """
        for dir_path in itertools.chain(dt_rpath, self.ld_library_path, dt_runpath, self.ld_so_paths):
            candidate_path = os.path.join(dir_path, lib_name)
            if os.path.exists(candidate_path):
                return candidate_path
        return None

    def ldd_shallow(self, path: str) -> List[Tuple[str, Optional[str]]]:
        """
        Perform ldd on the given binary and return a list of resolved library paths.
        """
        required_libs = []
        dt_rpath = []
        dt_runpath = []
        resolved_libs = []
        with open(path, "rb") as f:
            elffile = ELFFile(f)

            # Look for the dynamic section
            for section in elffile.iter_sections():
                if isinstance(section, DynamicSection):
                    # DT_NEEDED entries contain the names of required libraries
                    for tag in section.iter_tags():
                        if tag.entry.d_tag == "DT_NEEDED":
                            required_libs.append(tag.needed)
                        elif tag.entry.d_tag == "DT_RPATH":
                            dt_rpath.append(tag.rpath)
                        # RUNPATH (New style)
                        elif tag.entry.d_tag == "DT_RUNPATH":
                            dt_runpath.append(tag.runpath)
        dt_rpath = self._filter_ld_path(dt_rpath)
        dt_runpath = self._filter_ld_path(dt_runpath)
        for lib in required_libs:
            resolved_path = self._search_library_in_paths(lib, dt_rpath, dt_runpath)
            resolved_libs.append((lib, resolved_path))
        return resolved_libs

    def ldd_full(self, path: str) -> Dict[str, List[Tuple[str, Optional[str]]]]:
        """
        Perform full recursive ldd on the given binary and return a list of resolved library paths.
        """
        ldd_cache = {}
        # Perform a depth-first search to resolve all dependencies
        ld_dequeue = queue.SimpleQueue()
        ld_dequeue.put(path)
        for ld_preload in self.ld_preload:
            ld_dequeue.put(ld_preload)
        while not ld_dequeue.empty():
            current_path = ld_dequeue.get()
            if current_path in ldd_cache:
                continue
            resolved_libs = self.ldd_shallow(current_path)
            ldd_cache[current_path] = resolved_libs
            for lib_name, lib_path in resolved_libs:
                if lib_path is not None and lib_path not in ldd_cache:
                    ld_dequeue.put(lib_path)
        return ldd_cache

    def ldd_full_flattened(self, path: str) -> List[Tuple[str, Optional[str]]]:
        """
        Perform full recursive ldd on the given binary and return a flattened list of resolved library paths.
        """
        ldd_cache = self.ldd_full(path)
        flattened_list = []
        for libs in ldd_cache.values():
            flattened_list.extend(libs)
        return list(set(flattened_list))


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
        self.WITH_NCBI_NGS = False

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
        if self.WITH_NCBI_NGS:
            cmake_flags.append("-DWITH_NCBI_NGS=ON")
        if self.FIND_RANDOM_MKL_THROUGH_PKGCONF:
            cmake_flags.append(f"-DFIND_RANDOM_MKL_THROUGH_PKGCONF={self.FIND_RANDOM_MKL_THROUGH_PKGCONF}")
        return cmake_flags

    def generate_link_patterns(self) -> List[List[re.Pattern]]:
        """
        Generate link patterns for this build configuration.
        The 1st level of list is AND, the 2nd level is OR.

        :return:
        """
        retl = []
        if self.USE_MALLOC == "JEMALLOC":
            retl.append([re.compile("libjemalloc\\.so")])
        elif self.USE_MALLOC == "MIMALLOC":
            retl.append([re.compile("libmimalloc\\.so")])
        elif self.USE_MALLOC == "TCMALLOC":
            retl.append([re.compile("libtcmalloc\\.so")])
        elif self.USE_MALLOC == "TCMALLOC_MINIMAL":
            retl.append([re.compile("libtcmalloc_minimal\\.so")])

        if self.WITH_NCBI_NGS:
            retl.append([re.compile("libncbi-ngs\\.so")])

        if self.USE_HTSLIB is not None:
            retl.append([re.compile(f"lib{self.USE_HTSLIB}\\.so")])
        else:
            retl.append([re.compile("liblabw_slim_htslib\\.so")])

        if self.USE_LIBFMT is not None:
            retl.append([re.compile(f"lib{self.USE_LIBFMT}\\.so")])
        else:
            retl.append([re.compile("libslim_libfmt\\.so")])

        if self.USE_RANDOM_GENERATOR == "ONEMKL":
            retl.append(
                [
                    re.compile("libmkl_rt\\.so"),
                    re.compile("libmkl_intel_lp64\\.so"),
                    re.compile("libmkl_intel_ilp64\\.so"),
                ]
            )

        retl.append([re.compile("libboost_program_options\\.so")])
        retl.append([re.compile("libboost_log\\.so")])
        retl.append([re.compile("libboost_filesystem\\.so")])
        retl.append([re.compile("libboost_timer\\.so")])
        retl.append([re.compile("libz\\.so")])
        retl.append([re.compile("libslim_libceu\\.so")])
        return retl

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
        new_config.WITH_NCBI_NGS = self.WITH_NCBI_NGS
        return new_config


def do_build(config: BuildConfig, this_job_id: int) -> None:
    cmake_opts = config.generate_cmake_opts()
    os.makedirs(os.path.join(LOG_DIR, str(this_job_id)))
    build_dir = os.path.join(LOG_DIR, str(this_job_id), f"build")
    install_dir = os.path.join(LOG_DIR, str(this_job_id), f"install")
    art_modern_path = os.path.join(install_dir, "bin", "art_modern" + ("-mpi" if WITH_MPI else ""))
    apb_path = os.path.join(install_dir, "bin", "art_profile_builder" + ("-mpi" if WITH_MPI else ""))
    do_build_lh = logging.getLogger(f"DO_BUILD {this_job_id}/{num_total_jobs}")
    do_build_lh.handlers = []
    log_path = os.path.join(LOG_DIR, str(this_job_id), "main.log")
    do_build_lh.addHandler(logging.FileHandler(log_path, "w", "utf-8"))
    do_build_lh.setLevel(logging.INFO)
    for handler_ in do_build_lh.handlers:
        handler_.setLevel(logging.INFO)
        handler_.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    do_build_lh.info("Build options: %s", {" ".join(cmake_opts)})
    do_build_lh.info("B: %s, I: %s", build_dir, install_dir)
    do_build_lh.info("Build log for job %s", this_job_id)
    main_lh.info("%s/%s: Build options: %s", this_job_id, num_total_jobs, " ".join(cmake_opts))
    if not DRY_RUN:
        os.makedirs(build_dir)
        os.makedirs(install_dir)

    def failed():
        with IO_MUTEX:
            FAILED_ID.append(this_job_id)
        main_lh.error("%s/%s: FAILED", this_job_id, num_total_jobs)
        do_build_lh.error("FAILED")

    def run_wrapper(
        step_name_: str,
        cmdline: List[str],
        *args,
        additional_env: Optional[Dict[str, str]] = None,
        timeout: Optional[float] = None,
        **kwargs,
    ) -> bool:
        do_build_lh.info("%s START", step_name_)
        do_build_lh.info("%s CMDLINE: %s", step_name_, " ".join(cmdline))
        env = dict(os.environ)
        if additional_env:
            env.update(additional_env)
        with open(
            os.path.join(LOG_DIR, str(this_job_id), f"{step_name_.lower()}.env.sh"), "wt", encoding="utf-8"
        ) as env_file:
            for k, v in env.items():
                env_file.write(f"{k}='{v}'\n")
        if DRY_RUN:
            do_build_lh.info("%s DRY RUN", step_name_)
        else:
            with open(
                os.path.join(LOG_DIR, str(this_job_id), f"{step_name_.lower()}.log"), "wb", buffering=0
            ) as log_file:
                try:
                    proc = subprocess.Popen(cmdline, *args, env=env, stdout=log_file, stderr=log_file, **kwargs)
                except subprocess.SubprocessError as e:
                    do_build_lh.error("%s SUBPROCESS ERROR: %s", step_name_, str(e))
                    return False
                try:
                    if timeout is not None:
                        proc.wait(timeout=timeout)
                    else:
                        proc.wait()
                except subprocess.TimeoutExpired as e:
                    proc.terminate()
                    try:
                        proc.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    do_build_lh.error("%s TIMEOUT: %s", step_name_, str(e))
                    return False
                if proc.returncode != 0:
                    do_build_lh.error("%s FAILED with return code %s", step_name_, proc.returncode)
                    return False
        do_build_lh.info("%s DONE", step_name_)
        return True

    def run_ldd() -> bool:
        step_name_ = "LDD"
        do_build_lh.info("%s START", step_name_)
        ldd_lh = logging.getLogger(f"LDD {this_job_id}/{num_total_jobs}")
        ldd_lh.handlers = []
        ldd_lh.setLevel(logging.INFO)
        ldd_lh.addHandler(
            logging.FileHandler(os.path.join(LOG_DIR, str(this_job_id), f"{step_name_.lower()}.log"), "w", "utf-8")
        )
        for ldd_lh_handler_ in ldd_lh.handlers:
            ldd_lh_handler_.setLevel(logging.INFO)
            ldd_lh_handler_.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
        link_patterns = config.generate_link_patterns()
        if DRY_RUN:
            do_build_lh.info("%s DRY RUN", step_name_)
            ldd_lh.info(
                "Patterns: %s", "&".join("(" + "|".join(map(lambda x: x.pattern, lp)) + ")" for lp in link_patterns)
            )
        else:
            libs_flattened = set(
                itertools.chain(pls.ldd_full_flattened(art_modern_path), pls.ldd_full_flattened(apb_path))
            )
            for lib, path in libs_flattened:
                ldd_lh.info("LIB %s => %s", lib, path)
            missing_libs = [lib for lib, path in libs_flattened if path is None]
            if missing_libs:
                ldd_lh.info("MISSING LIBS: %s", ", ".join(missing_libs))
                return False
            link_patterns = config.generate_link_patterns()
            for link_pattern_group in link_patterns:
                ldd_pattern_str = " | ".join(p.pattern for p in link_pattern_group)
                found = False
                for lib, path in libs_flattened:
                    if found:
                        break
                    for p in link_pattern_group:
                        if p.search(path):
                            ldd_lh.info("FOUND: %s matches %s", lib, p.pattern)
                            found = True
                            break
                if not found:
                    ldd_lh.info("FAILED: Nothing match %s", ldd_pattern_str)
                    return False
        ldd_lh.info("DONE")
        do_build_lh.info("%s DONE", step_name_)
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
    ):
        failed()
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
        timeout=360,
    ):
        failed()
        return
    if not run_wrapper(
        "CTEST",
        [CTEST, "--output-on-failure"],
        cwd=build_dir,
        timeout=60,
    ):
        failed()
        return

    if not run_wrapper(
        "INSTALL",
        [CMAKE, "--install", build_dir],
        timeout=60,
    ):
        failed()
        return

    if not run_ldd():
        failed()
        return

    if not run_wrapper(
        "TESTSMALL",
        [BASH, os.path.join(SHDIR, "test-small.sh")],
        additional_env={
            **os.environ,
            "ART_MODERN_PATH": art_modern_path,
            "APB_PATH": apb_path,
            "HELP_VERSION_ONLY": "1" if config.HELP_VERSION_ONLY else "0",
            "MPIEXEC": MPIEXEC if WITH_MPI else "",
            "SAMTOOLS_THREADS": "4",
            "MPI_PARALLEL": "4",
            "PARALLEL": "2",
            "NO_FASTQC": "1",
            "OUT_DIR": os.path.join(LOG_DIR, str(this_job_id), f"testmall.out.d"),
            "SET_X": "1",
            "WITH_NCBI_NGS": "1" if config.WITH_NCBI_NGS else "0",
        },
        timeout=1200,
    ):
        failed()
        return
    if not DRY_RUN:
        shutil.rmtree(build_dir)
        shutil.rmtree(install_dir)
        shutil.rmtree(os.path.join(LOG_DIR, str(this_job_id)))
    with IO_MUTEX:
        SUCCESS_ID.append(this_job_id)
    do_build_lh.info("SUCCESS")
    main_lh.info("%s/%s: SUCCESS", this_job_id, num_total_jobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test build script re-implementation.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a dry run without executing commands.",
    )
    parser.add_argument(
        "--ld-run",
        action="store_true",
        help="Stop after parsing ld.so.conf and LD_LIBRARY_PATH.",
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
    parser.add_argument(
        "--no-clear-on-success",
        action="store_true",
        help="Do not automatically clear log dir on success.",
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
    main_lh = logging.getLogger("Main")
    main_lh.setLevel(logging.INFO)
    main_lh.addHandler(logging.StreamHandler(sys.stderr))
    main_lh.addHandler(logging.FileHandler(os.path.join(LOG_DIR, "main.log")))
    for handler in main_lh.handlers:
        handler.setLevel(logging.INFO)
        handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    main_lh.info("Log dir: %s", LOG_DIR)
    pls = PyLddSearch()
    if _args.ld_run:
        exit(0)

    RANDOM_GENERATORS = ["STL", "PCG", "BOOST"]

    print(f"Probing for MKL through CMake...", end="")
    if probe_using_cmake_project("test-mkl.cmake"):
        RANDOM_GENERATORS.append("ONEMKL")
        print("SUCCESS")
    else:
        print("FAIL")

    print(f"Probing for concurrent queue...", end="")
    concurrent_queue_path = probe_concurrent_queue()
    if concurrent_queue_path:
        print(concurrent_queue_path)
    else:
        print("FAIL")

    MALLOCS = ["NOP"]
    print(f"Probing for jemalloc...", end="")
    if probe_pkg_config("jemalloc"):
        MALLOCS.append("JEMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print(f"Probing for mimalloc...", end="")
    if probe_using_cmake_project("test-mimalloc.cmake"):
        MALLOCS.append("MIMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print(f"Probing for tcmalloc...", end="")
    if probe_pkg_config("libtcmalloc"):
        MALLOCS.append("TCMALLOC")
        print("SUCCESS")
    else:
        print("FAIL")
    print(f"Probing for tcmalloc_minimal...", end="")
    if probe_pkg_config("libtcmalloc_minimal"):
        MALLOCS.append("TCMALLOC_MINIMAL")
        print("SUCCESS")
    else:
        print("FAIL")

    print(f"Probing for HTSLib...", end="")
    is_htslib_exist = probe_pkg_config("htslib")
    if is_htslib_exist:
        print("SUCCESS")
    else:
        print("FAIL")

    print(f"Probing for " + "{fmt}" + "...", end="")
    is_libfmt_exist = probe_pkg_config("fmt")
    if is_libfmt_exist:
        print("SUCCESS")
    else:
        print("FAIL")

    print(f"Probing for <pcg_random.hpp>...", end="")
    pcg_random_hpp_exist_path = probe_pcg_random_hpp()
    if pcg_random_hpp_exist_path:
        RANDOM_GENERATORS.append("SYSTEM_PCG")
        print(pcg_random_hpp_exist_path)
    else:
        print("FAIL")

    print(f"Probing for MKL found through pkg-config...", end="")
    mkl_pc = peobe_mkl_using_pkgconf()
    if mkl_pc:
        print(";".join(mkl_pc))
    else:
        print("FAIL")

    print("Probing for NCBI NGS libraries...", end="")
    have_ncbi_ngs = probe_ncbi_ngs()
    if have_ncbi_ngs:
        print("SUCCESS")
    else:
        print("FAIL")

    cmake_build_types = ["Debug", "Release", "RelWithDebInfo"]
    max_workers = min(5, os.cpu_count() // 4)
    bcs = []

    for cmake_build_type in cmake_build_types:
        bc = BuildConfig(_old_cmake_flags)

        bc.HELP_VERSION_ONLY = not _args.small
        bc.CMAKE_BUILD_TYPE = cmake_build_type

        for use_thread_parallel in ["BS", "NOP"]:
            bc.USE_THREAD_PARALLEL = use_thread_parallel
            bcs.append(bc.copy())
        bc.USE_THREAD_PARALLEL = "ASIO"

        for use_malloc in MALLOCS:
            bc.USE_MALLOC = use_malloc
            bcs.append(bc.copy())
        bc.USE_MALLOC = "AUTO"

        if is_libfmt_exist:
            bc.USE_LIBFMT = "fmt"
            bcs.append(bc.copy())
            bc.USE_LIBFMT = None

        if is_htslib_exist:
            bc.USE_HTSLIB = "hts"
            bcs.append(bc.copy())
            bc.USE_HTSLIB = None

        if concurrent_queue_path is not None:
            bc.USE_CONCURRENT_QUEUE = concurrent_queue_path
            bcs.append(bc.copy())
            bc.USE_CONCURRENT_QUEUE = None

        bc.AM_NO_Q_REVERSE = True
        bcs.append(bc.copy())
        bc.AM_NO_Q_REVERSE = False

        for use_random_generator in RANDOM_GENERATORS:
            bc.USE_RANDOM_GENERATOR = use_random_generator
            bcs.append(bc.copy())
        if mkl_pc:
            bc.USE_RANDOM_GENERATOR = "ONEMKL"
            for pc in mkl_pc:
                bc.FIND_RANDOM_MKL_THROUGH_PKGCONF = pc
                bcs.append(bc.copy())
            bc.FIND_RANDOM_MKL_THROUGH_PKGCONF = None

        bc.USE_RANDOM_GENERATOR = "PCG"

        if have_ncbi_ngs:
            bc.WITH_NCBI_NGS = True
            bcs.append(bc.copy())
            bc.WITH_NCBI_NGS = False
    num_total_jobs = len(bcs)
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for job_id, bc in enumerate(bcs, 1):
            futures.append(executor.submit(do_build, bc, job_id))
        for future in futures:
            try:
                _ = future.result()
            except concurrent.futures.CancelledError:
                # Ignore
                pass
            except Exception as e:
                main_lh.error("Exception inside TPE: %s: %s", str(e.__class__.__name__), str(e))

                executor.shutdown(cancel_futures=True)
    SUCCESS_ID = sorted(SUCCESS_ID)
    FAILED_ID = sorted(FAILED_ID)
    main_lh.info("Successful builds: %s", ", ".join(map(str, SUCCESS_ID)))
    main_lh.info("Failed builds: %s", ", ".join(map(str, FAILED_ID)))
    if not FAILED_ID and not _args.no_clear_on_success:
        shutil.rmtree(LOG_DIR)
    exit(1 if FAILED_ID else 0)
