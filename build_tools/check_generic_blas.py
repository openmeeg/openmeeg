"""Check a BLA_IMPLEMENTATION=Generic build of OpenMEEG.

Usage::

    python check_generic_blas.py linkage <path/to/OpenMEEGMaths library>
    python check_generic_blas.py runtime <netlib|openblas|mkl|accelerate>

``linkage`` asserts that the library links only to the generic
BLAS/CBLAS/LAPACK/LAPACKE libraries and not to a specific implementation, so
that the implementation can be swapped at install time. ``runtime`` asserts
which implementation actually gets loaded alongside OpenMEEG.
"""

import re
import subprocess
import sys
from pathlib import Path

GENERIC = ("blas", "cblas", "lapack", "lapacke")
SPECIFIC = ("openblas", "mkl", "blis", "flexiblas", "accelerate", "veclib")


def _dependencies(lib):
    if sys.platform == "win32":
        import pefile

        pe = pefile.PE(lib, fast_load=True)
        pe.parse_data_directories(
            directories=[pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"]]
        )
        return [entry.dll.decode() for entry in pe.DIRECTORY_ENTRY_IMPORT]
    if sys.platform == "darwin":
        out = subprocess.check_output(["otool", "-L", lib], text=True)
        return [line.split()[0] for line in out.splitlines()[1:]]
    out = subprocess.check_output(["readelf", "-d", lib], text=True)
    return re.findall(r"\(NEEDED\).*\[(.+)\]", out)


def _stem(dep):
    # libcblas.so.3 / @rpath/libcblas.3.dylib / libcblas.dll -> cblas
    name = Path(dep).name.lower()
    name = name.removeprefix("lib")
    return name.split(".")[0]


def check_linkage(lib):
    """Assert that lib links to the generic libraries only."""
    deps = _dependencies(lib)
    print(f"{lib} depends on:\n  " + "\n  ".join(deps))
    stems = {_stem(dep) for dep in deps}
    missing = set(GENERIC) - stems
    assert not missing, f"Not linked to the generic {sorted(missing)}"
    specific = [dep for dep in deps if any(s in dep.lower() for s in SPECIFIC)]
    assert not specific, f"Linked to a specific implementation: {specific}"


def _accelerate_loaded():
    # Accelerate's BLAS/LAPACK live in vecLib, which threadpoolctl does not
    # report, so look at what dyld has loaded into this process instead
    import ctypes

    libc = ctypes.CDLL(None)
    libc._dyld_image_count.restype = ctypes.c_uint32
    libc._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    libc._dyld_get_image_name.restype = ctypes.c_char_p
    names = [
        libc._dyld_get_image_name(ii).decode() for ii in range(libc._dyld_image_count())
    ]
    return any("/vecLib.framework/" in name for name in names)


def check_runtime(expected):
    """Assert that the expected implementation gets loaded with OpenMEEG."""
    import threadpoolctl

    import openmeeg  # noqa: F401  -- importing loads the extension and its BLAS

    infos = [
        info for info in threadpoolctl.threadpool_info() if info["user_api"] == "blas"
    ]
    for info in infos:
        print(f"{info['internal_api']}: {info['filepath']}")
    got = {info["internal_api"] for info in infos}
    # threadpoolctl recognizes neither the reference implementation nor Accelerate
    want = set() if expected in ("netlib", "accelerate") else {expected}
    assert got == want, f"Expected BLAS implementation {want}, got {got}"
    if sys.platform == "darwin":
        accelerate = _accelerate_loaded()
        print(f"Accelerate loaded: {accelerate}")
        assert accelerate == (expected == "accelerate")


if __name__ == "__main__":
    mode, arg = sys.argv[1:]
    {"linkage": check_linkage, "runtime": check_runtime}[mode](arg)
