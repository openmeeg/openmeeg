"""Adapted from SciPy."""

import os

# arm64 macOS needed OPENBLAS_THREAD_TIMEOUT=1 here (scipy#15050) until those
# wheels moved to Accelerate and stopped shipping an OpenBLAS to clash with
if os.name == "nt":
    import glob
    from ctypes import WinDLL

    libs_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    try:
        owd = os.getcwd()
        os.chdir(libs_path)
        for filename in glob.glob(os.path.join(libs_path, ".libs", "*dll")):
            WinDLL(os.path.abspath(filename))
    finally:
        os.chdir(owd)
