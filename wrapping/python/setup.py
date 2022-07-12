#!/usr/bin/env python

# Copyright (C) 2011-2020 Alexandre Gramfort
# <alexandre.gramfort@inria.fr>
# This assumes/requires that SWIG has already been run so the .so exists!

from pathlib import Path
import os
from setuptools import setup, Extension  # noqa

root = Path(__file__).parent

version = None
with open((root / 'openmeeg' / '_version.py'), 'r') as fid:
    for line in (line.strip() for line in fid):
        if line.startswith('__version__'):
            version = line.split('=')[1].strip().strip('\'')
            break
if version is None:
    raise RuntimeError('Could not determine version')


DISTNAME = 'openmeeg'
DESCRIPTION = 'Forward problems solver in the field of EEG and MEG.'
MAINTAINER = 'Alexandre Gramfort'
MAINTAINER_EMAIL = 'alexandre.gramfort@inria.fr'
URL = 'https://openmeeg.github.io/'
LICENSE = 'CECILL-B'
DOWNLOAD_URL = 'http://github.com/openmeeg/openmeeg'
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


"""  # in case someday we do a proper SWIG build...
from ctypes.util import find_library
import glob
import platform
from shutil import copyfile
import numpy as np
from distutils.command.build import build


def find_openmeeg_include():
    # If a path is provded by the environment, expect the GLPK header there.
    header_path = os.environ.get('OPENMEEG_HEADER_PATH', None)
    if not header_path:
        header_path = find_library('OpenMEEG')
        if header_path is not None:
            header_path = Path(header_path).parent.parent / 'include'
    else:
        header_path = Path(header_path)

    if header_path is None or not header_path.is_dir():
        extra = f' {header_path}' if header_path is not None else ''
        raise RuntimeError(f'Could not find OpenMEEG header directory{extra}, '
                           'consider setting OPENMEEG_HEADER_PATH')
    header_path = header_path.resolve()
    header_file = header_path / 'OpenMEEG' / 'OpenMEEGMathsConfig.h'
    if not header_file.is_file():
        raise RuntimeError(
            f'Could not find expected OpenMEEG header {header_file} in '
            f'directory {header_path}')
    return header_path
"""

if __name__ == "__main__":
    import numpy as np
    manifest = (root / 'MANIFEST')
    if manifest.is_file():
        os.remove(manifest)

    with open('README.rst', 'r') as fid:
        long_description = fid.read()

    # SWIG
    ext_modules = []
    if os.getenv('OPENMEEG_USE_SWIG', '0').lower() in ('1', 'true'):
        include_dirs = [np.get_include()]
        swig_opts = ['-c++', '-v', '-Werror']
        library_dirs = []
        openmeeg_include = os.getenv('OPENMEEG_INCLUDE')
        if openmeeg_include is not None:
            openmeeg_include = Path(openmeeg_include)
            assert openmeeg_include.is_dir(), openmeeg_include
            include_dirs.append(str(openmeeg_include))
            swig_opts.append(f'-I{openmeeg_include}')
        openmeeg_lib = os.getenv('OPENMEEG_LIB')
        if openmeeg_lib is not None:
            openmeeg_lib = Path(openmeeg_lib)
            assert openmeeg_lib.is_dir(), openmeeg_lib
            library_dirs.append(openmeeg_lib)
        swig_openmeeg = Extension(
            "openmeeg._openmeeg",
            ["openmeeg/openmeeg.i"],
            libraries=['OpenMEEG'],
            swig_opts=swig_opts,
            extra_compile_args=['-v', '-std=c++17'],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
        )
        ext_modules.append(swig_openmeeg)

    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          include_package_data=True,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          long_description=long_description,
          long_description_content_type='text/x-rst',
          zip_safe=False,  # the package can run out of an .egg file
          classifiers=['Intended Audience :: Science/Research',
                       'Intended Audience :: Developers',
                       'License :: OSI Approved',
                       'Programming Language :: Python',
                       'Topic :: Software Development',
                       'Topic :: Scientific/Engineering',
                       'Operating System :: Microsoft :: Windows',
                       'Operating System :: POSIX',
                       'Operating System :: Unix',
                       'Operating System :: MacOS',
                       'Programming Language :: Python :: 3',
                       ],
          keywords='neuroscience neuroimaging MEG EEG ECoG sEEG iEEG brain',
          project_urls={
              'Documentation': 'https://openmeeg.github.io',
              'Source': 'https://github.com/openmeeg/openmeeg',
              'Tracker': 'https://github.com/openmeeg/openmeeg/issues',
          },
          platforms='any',
          python_requires='>=3.7',
          install_requires=["numpy"],
          packages=["openmeeg", "openmeeg.tests"],
          #cmdclass={  # TODO: This breaks macOS for some reason!
          #    'bdist_wheel': bdist_wheel,
          #},
          ext_modules=ext_modules,
          )
