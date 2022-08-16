"""Adapted from SciPy."""

import os
import sys
import platform

if os.name == "nt":
    from ctypes import WinDLL
    import glob

    libs_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    try:
        owd = os.getcwd()
        os.chdir(libs_path)
        for filename in glob.glob(os.path.join(libs_path, ".libs", "*dll")):
            WinDLL(os.path.abspath(filename))
    finally:
        os.chdir(owd)
elif sys.platform == "darwin" and platform.machine() == "arm64":
    # On arm64 macOS the OpenBLAS runtimes of NumPy and SciPy don't seem to work
    # well together unless this timeout limit is set - it results in very slow
    # performance for some linalg functionality.
    # See https://github.com/scipy/scipy/issues/15050 for details.
    os.environ["OPENBLAS_THREAD_TIMEOUT"] = "1"
