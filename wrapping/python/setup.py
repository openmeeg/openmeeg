#!/usr/bin/env python

# Copyright (C) 2011-2020 Alexandre Gramfort
# <alexandre.gramfort@inria.fr>

from pathlib import Path
import os
import platform
import sys
import sysconfig

from setuptools import setup, Extension  # noqa
from setuptools.command import build_py
from setuptools.dist import Distribution


# Oldest interpreter the stable ABI wheel targets, as a version hex and as the
# matching wheel tag. The free-threaded build has no limited API at all --
# Python.h #errors out if Py_LIMITED_API is defined alongside Py_GIL_DISABLED,
# and setuptools raises on py_limited_api there too (CPython gh-111506).
PY_LIMITED_API = "0x030A0000"
PY_LIMITED_API_TAG = "cp310"
abi3 = (
    platform.python_implementation() == "CPython"
    and not sysconfig.get_config_var("Py_GIL_DISABLED")
)


# CMake-driven builds hand us a prebuilt extension rather than an ext_module, so
# setuptools would otherwise call the wheel pure -- which loses the platform tag
# (nicer for app building) and suppresses the abi3 tag entirely.
class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


# Subclass the build command so that build_ext is called before build_py
class BuildExtFirst(build_py.build_py):
    def run(self):
        self.run_command("build_ext")
        super().run()


if __name__ == "__main__":
    import numpy as np

    # SWIG
    cmdclass = dict(build_py=BuildExtFirst)
    ext_modules = []
    if os.getenv("OPENMEEG_USE_SWIG", "0").lower() in ("1", "true"):
        include_dirs = [np.get_include()]
        swig_opts = [
            "-c++",
            "-v",
            "-module",
            "_openmeeg_wrapper",
            "-interface",
            "_openmeeg",
            # "-O" without its other two flags: -fastproxy emits PyCFunctionObject and
            # PyInstanceMethod_New, neither in the limited API, and -fvirtual measured
            # ~6.5% slower here (it routes the 8 LinOp-derived size()/info() through a
            # base-class ConvertPtr). -fastdispatch is purely an optimization now that
            # the typecheck typemaps in _openmeeg.i make overload resolution work
            # without it, so it can be dropped if upstream ever stops supporting it
            # alongside abi3 (swig/swig#3089).
            "-fastdispatch",
            # Declares Py_MOD_GIL_NOT_USED (SWIG >= 4.4, hence the build-system pin),
            # which the i/o layer earned in gh-866. Without it CPython silently
            # re-enables the GIL on free-threaded builds.
            "-nogil",
            # Someday we could look at other options like:
            # "-extranative",  # Return extra native wrappers for C++ std containers wherever possible
            # "-castmode",  # Enable the casting mode, which allows implicit cast between types in Python
            # "-flatstaticmethod",  # Generate additional flattened Python methods for C++ static methods
            # "-olddefs",  # Keep the old method definitions when using -fastproxy
            # TODO someday we should add:
            # "-Werror",
        ]
        library_dirs = []
        # Read from the installed config header rather than re-deriving the
        # per-platform rules here, so this tracks whatever cmake actually built
        use_accelerate = False
        openmeeg_include = os.getenv("OPENMEEG_INCLUDE")
        if openmeeg_include is not None:
            openmeeg_include = Path(openmeeg_include).resolve(strict=True)
            print(f"Adding OpenMEEG include directory: {openmeeg_include}")
            include_dirs.append(str(openmeeg_include))
            swig_opts.append(f"-I{openmeeg_include}")
            config_h = openmeeg_include / "OpenMEEGConfigure.h"
            if config_h.is_file():
                use_accelerate = "#define USE_ACCELERATE" in config_h.read_text()
            print(f"Building against Accelerate: {use_accelerate}")
        msvc = os.getenv("SWIG_FLAGS", "") == "msvc"
        # Empty means "not applicable" (Accelerate builds have no OpenBLAS tree),
        # since cibuildwheel's per-arch overrides can only override, not unset.
        openblas_include = os.getenv("OPENBLAS_INCLUDE") or None
        if openblas_include is not None:
            openblas_include = Path(openblas_include).resolve(strict=True)
            print(f"Adding OpenBLAS include directory: {openblas_include}")
            include_dirs.append(str(openblas_include))
            swig_opts.append(f"-I{openblas_include}")
        openmeeg_lib = os.getenv("OPENMEEG_LIB")
        if openmeeg_lib is not None:
            openmeeg_lib = Path(openmeeg_lib).resolve(strict=True)
            print(f"Adding OpenMEEG library directory: {openmeeg_lib}")
            library_dirs.append(str(openmeeg_lib))
        extra_compile_opts, extra_link_opts = [], []
        if msvc:
            extra_compile_opts.extend(["/std:c++17"])
            extra_link_opts.append("OpenMEEGMaths.lib")
        else:
            extra_compile_opts.extend(["-v", "-std=c++17"])
        if sys.platform == "darwin":
            if use_accelerate:
                version_min = "13.3"  # first release with Accelerate's new BLAS/LAPACK
            elif "arm64" in platform.platform():
                version_min = "11"
            else:
                version_min = "10.15"
            extra_compile_opts.extend([f"-mmacosx-version-min={version_min}"])
            if use_accelerate:
                # The inline BLAS calls in matrix.h land in this extension too, so
                # don't leave them to flat-namespace lookup via libOpenMEEG
                extra_link_opts.append("-framework")
                extra_link_opts.append("Accelerate")
        # An example cmake command that works is:
        #   C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Tools\MSVC\14.16.27023\bin\HostX64\x64\link.exe
        #     /ERRORREPORT:QUEUE /INCREMENTAL:NO /NOLOGO
        #     /OUT:"D:\a\openmeeg\openmeeg\build\wrapping\python\openmeeg\_openmeeg.pyd"
        #     C:\Miniconda\envs\test\libs\python310.lib ..\..\OpenMEEG\Release\OpenMEEG.lib ..\..\OpenMEEGMaths\Release\OpenMEEGMaths.lib ..\..\..\openblas\64\lib\openblas.lib "..\..\vcpkg_installed\x64-windows-release\lib\libmatio.lib" "..\..\vcpkg_installed\x64-windows-release\lib\hdf5.lib" kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib
        #     /MANIFEST /MANIFESTUAC:"level='asInvoker' uiAccess='false'" /manifest:embed /PDB:"D:/a/openmeeg/openmeeg/build/wrapping/python/openmeeg/_openmeeg.pdb"
        #     /SUBSYSTEM:CONSOLE /TLBID:1 /DYNAMICBASE /NXCOMPAT
        #     /IMPLIB:"D:/a/openmeeg/openmeeg/build/wrapping/python/Release/openmeeg.lib" /MACHINE:X64 /DLL openmeeg.dir\Release\openmeegPYTHON_wrap.obj
        # An example setuptools command that fails is:
        #   C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Tools\MSVC\14.29.30133\bin\HostX86\x64\cl.exe
        #     -IC:\Miniconda\envs\test\lib\site-packages\numpy\core\include -ID:\a\openmeeg\openmeeg\install\include\OpenMEEG -ID:\a\openmeeg\openmeeg\openblas\64\include -IC:\Miniconda\envs\test\include -IC:\Miniconda\envs\test\Include "-IC:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Tools\MSVC\14.29.30133\ATLMFC\include" "-IC:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Tools\MSVC\14.29.30133\include" "-IC:\Program Files (x86)\Windows Kits\NETFXSDK\4.8\include\um" "-IC:\Program Files (x86)\Windows Kits\10\include\10.0.22000.0\ucrt" "-IC:\Program Files (x86)\Windows Kits\10\include\10.0.22000.0\shared" "-IC:\Program Files (x86)\Windows Kits\10\include\10.0.22000.0\um" "-IC:\Program Files (x86)\Windows Kits\10\include\10.0.22000.0\winrt" "-IC:\Program Files (x86)\Windows Kits\10\include\10.0.22000.0\cppwinrt"
        #     /c /nologo /O2 /W3 /GL /DNDEBUG /MD
        #     /EHsc /Tpopenmeeg/openmeeg_wrap.cpp /Fobuild\temp.win-amd64-cpython-310\Release\openmeeg/openmeeg_wrap.obj

        define_macros = []
        abi3_kwargs = dict()
        if abi3:
            define_macros += [
                ("Py_LIMITED_API", PY_LIMITED_API),
            ]
        swig_openmeeg = Extension(
            "openmeeg._openmeeg",
            sources=["openmeeg/_openmeeg.i"],
            libraries=["OpenMEEG"],
            swig_opts=swig_opts,
            define_macros=define_macros,
            extra_compile_args=extra_compile_opts,
            extra_link_args=extra_link_opts,
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            py_limited_api=abi3,
        )
        ext_modules.append(swig_openmeeg)

    setup(
        cmdclass=cmdclass,
        distclass=BinaryDistribution,
        ext_modules=ext_modules,
        options=(
            {"bdist_wheel": {"py_limited_api": PY_LIMITED_API_TAG}} if abi3 else {}
        ),
    )
