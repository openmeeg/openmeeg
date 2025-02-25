#!/usr/bin/env python

# Copyright (C) 2011-2020 Alexandre Gramfort
# <alexandre.gramfort@inria.fr>

from pathlib import Path
import os
import platform
import sys

from setuptools import setup, Extension  # noqa
from setuptools.command import build_py
from wheel.bdist_wheel import bdist_wheel


abi3 = (platform.python_implementation() == "CPython")
abi3 = abi3 and not os.getenv("OPENMEEG_NO_ABI3", "0").lower() in ("1", "true")


# Adapted from Apache-2.0 licensed code at:
# https://github.com/joerick/python-abi3-package-sample/blob/main/setup.py

class bdist_wheel_abi3(bdist_wheel):
    def get_tag(self):
        python, abi, plat = super().get_tag()
        if abi3:
            python, abi = "cp310", "abi3"
        return python, abi, plat

    # In cases where we don't SWIG, we still want to mark the wheel as impure
    # (to make things nicer for app building)
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False


# Subclass the build command so that build_ext is called before build_py
class BuildExtFirst(build_py.build_py):
    def run(self):
        self.run_command("build_ext")
        super().run()


if __name__ == "__main__":
    import numpy as np

    # SWIG
    cmdclass = dict(build_py=BuildExtFirst, bdist_wheel=bdist_wheel_abi3)
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
            # -O is problematic for abi3 because it enables "-fastproxy" which leads to
            # linking errors. However, "-fastdispatch" is needed to create our
            # typemaps, which seems like a bug (?) but doesn't create any ABI3 compat
            # issues.
            "-fastdispatch",
            # Someday we could look at other options like:
            # "-extranative",  # Return extra native wrappers for C++ std containers wherever possible
            # "-castmode",  # Enable the casting mode, which allows implicit cast between types in Python
            # "-flatstaticmethod",  # Generate additional flattened Python methods for C++ static methods
            # "-olddefs",  # Keep the old method definitions when using -fastproxy
            # TODO someday we should add:
            # "-Werror",
        ]
        library_dirs = []
        openmeeg_include = os.getenv("OPENMEEG_INCLUDE")
        if openmeeg_include is not None:
            openmeeg_include = Path(openmeeg_include).resolve(strict=True)
            include_dirs.append(str(openmeeg_include))
            swig_opts.append(f"-I{openmeeg_include}")
        msvc = os.getenv("SWIG_FLAGS", "") == "msvc"
        openblas_include = os.getenv("OPENBLAS_INCLUDE")
        if openblas_include is not None:
            openblas_include = Path(openblas_include).resolve(strict=True)
            include_dirs.append(str(openblas_include))
            swig_opts.append(f"-I{openblas_include}")
        openmeeg_lib = os.getenv("OPENMEEG_LIB")
        if openmeeg_lib is not None:
            openmeeg_lib = Path(openmeeg_lib).resolve(strict=True)
            library_dirs.append(str(openmeeg_lib))
        extra_compile_opts, extra_link_opts = [], []
        if msvc:
            extra_compile_opts.extend(["/std:c++17"])
            extra_link_opts.append("OpenMEEGMaths.lib")
        else:
            extra_compile_opts.extend(["-v", "-std=c++17"])
        if sys.platform == "darwin":
            version_min = "11" if "arm64" in platform.platform() else "10.15"
            extra_compile_opts.extend([f"-mmacosx-version-min={version_min}"])
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

        define_macros = [("SWIG_PYTHON_SILENT_MEMLEAK", None)]
        abi3_kwargs = dict()
        if abi3:
            define_macros += [
                ("Py_LIMITED_API", "0x030A0000"),  # 3.10
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
        ext_modules=ext_modules,
    )
