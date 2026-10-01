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

# OpenMEEG calls the C interfaces only, which in turn use blas and lapack. Those
# show up as direct dependencies too, except on Windows, where the linker drops
# DLLs that nothing calls into directly.
GENERIC = ("cblas", "lapacke")
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


def _accelerate_provides_cblas():
    # threadpoolctl does not report Accelerate. Merely finding it loaded proves
    # nothing either, since other libraries pull it in (e.g., HDF5's S3 support
    # via system frameworks), so ask which image actually defines the cblas_dgemm
    # that the generic libcblas OpenMEEG links to resolves to.
    import ctypes

    class DlInfo(ctypes.Structure):
        _fields_ = [
            ("dli_fname", ctypes.c_char_p),
            ("dli_fbase", ctypes.c_void_p),
            ("dli_sname", ctypes.c_char_p),
            ("dli_saddr", ctypes.c_void_p),
        ]

    cblas = ctypes.CDLL(str(Path(sys.prefix) / "lib" / "libcblas.3.dylib"))
    libc = ctypes.CDLL(None)
    libc.dladdr.argtypes = [ctypes.c_void_p, ctypes.POINTER(DlInfo)]
    info = DlInfo()
    addr = ctypes.cast(cblas.cblas_dgemm, ctypes.c_void_p)
    assert libc.dladdr(addr, ctypes.byref(info)), "dladdr failed"
    fname = info.dli_fname.decode()
    print(f"cblas_dgemm comes from {fname}")
    return "/vecLib.framework/" in fname


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
        assert _accelerate_provides_cblas() == (expected == "accelerate")


if __name__ == "__main__":
    mode, arg = sys.argv[1:]
    {"linkage": check_linkage, "runtime": check_runtime}[mode](arg)
