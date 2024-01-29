#!/usr/bin/env python

# Copyright (C) 2011-2020 Alexandre Gramfort
# <alexandre.gramfort@inria.fr>
# This assumes/requires that SWIG has already been run so the .so exists!

from pathlib import Path
import os
import platform
import sys

from setuptools import setup, Extension  # noqa
from setuptools.command import build_py

root = Path(__file__).parent

version = None
with open((root / "openmeeg" / "_version.py"), "r") as fid:
    for line in (line.strip() for line in fid):
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip('"')
            break
if version is None:
    raise RuntimeError("Could not determine version")


DISTNAME = "openmeeg"
DESCRIPTION = "Forward problems solver in the field of EEG and MEG."
MAINTAINER = "Alexandre Gramfort"
MAINTAINER_EMAIL = "alexandre.gramfort@inria.fr"
URL = "https://openmeeg.github.io/"
LICENSE = "CECILL-B"
DOWNLOAD_URL = "http://github.com/openmeeg/openmeeg"
VERSION = version

# Adapted from MIT-licensed
# https://github.com/Yelp/dumb-init/blob/48db0c0d0ecb4598d1a6400710445b85d67616bf/setup.py#L11-L27  # noqa
try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            # Mark us as not a pure python package
            self.root_is_pure = False

except ImportError:
    bdist_wheel = None  # noqa


# Subclass the build command so that build_ext is called before build_py
class BuildExtFirst(build_py.build_py):
    def run(self):
        self.run_command("build_ext")
        super().run()


if __name__ == "__main__":
    import numpy as np

    manifest = root / "MANIFEST"
    if manifest.is_file():
        os.remove(manifest)

    with open("README.rst", "r") as fid:
        long_description = fid.read()

    # SWIG
    cmdclass = dict(build_py=BuildExtFirst)
    ext_modules = []
    if os.getenv("OPENMEEG_USE_SWIG", "0").lower() in ("1", "true"):
        include_dirs = [np.get_include()]
        swig_opts = [
            "-c++",
            "-v",
            "-O",
            "-module",
            "_openmeeg_wrapper",
            "-interface",
            "_openmeeg",
            "-modern",
        ]  # TODO: , '-Werror']
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
            extra_link_opts.extend(["/std:c++17"])
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

        swig_openmeeg = Extension(
            "openmeeg._openmeeg",
            sources=["openmeeg/_openmeeg.i"],
            libraries=["OpenMEEG"],
            swig_opts=swig_opts,
            define_macros=[("SWIG_PYTHON_SILENT_MEMLEAK", None)],
            extra_compile_args=extra_compile_opts,
            include_dirs=include_dirs,
            library_dirs=library_dirs,
        )
        ext_modules.append(swig_openmeeg)
    else:  # built with -DENABLE_PYTHON=ON
        if sys.platform != "darwin" or os.getenv(
            "OPENMEEG_MACOS_WHEEL_PURE", "true"
        ).lower() in ("false", "0"):
            cmdclass["bdist_wheel"] = bdist_wheel

    setup(
        name=DISTNAME,
        maintainer=MAINTAINER,
        include_package_data=True,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        long_description=long_description,
        long_description_content_type="text/x-rst",
        zip_safe=False,  # the package can run out of an .egg file
        classifiers=[
            "Intended Audience :: Science/Research",
            "Intended Audience :: Developers",
            "License :: OSI Approved",
            "Programming Language :: Python",
            "Topic :: Software Development",
            "Topic :: Scientific/Engineering",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Operating System :: MacOS",
            "Programming Language :: Python :: 3",
        ],
        keywords="neuroscience neuroimaging MEG EEG ECoG sEEG iEEG brain",
        project_urls={
            "Documentation": "https://openmeeg.github.io",
            "Source": "https://github.com/openmeeg/openmeeg",
            "Tracker": "https://github.com/openmeeg/openmeeg/issues",
        },
        platforms="any",
        python_requires=">=3.9",
        install_requires=["numpy"],
        packages=["openmeeg", "openmeeg.tests"],
        cmdclass=cmdclass,
        ext_modules=ext_modules,
    )
