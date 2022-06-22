#!/usr/bin/env python

# Copyright (C) 2011-2020 Alexandre Gramfort
# <alexandre.gramfort@inria.fr>
# This assumes/requires that SWIG has already been run so the .so exists!

from ctypes.util import find_library
import glob
from pathlib import Path
import os
import platform  # noqa
from shutil import copyfile

import numpy as np
from setuptools import setup, Extension
from distutils.command.build import build

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


"""  # in case someday we do a proper SWIG build...
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
    manifest = (root / 'MANIFEST')
    if manifest.is_file():
        os.remove(manifest)

    with open('README.rst', 'r') as fid:
        long_description = fid.read()

    #openmeeg_root = root / '..' / '..' / 'install'
    #openmeeg_include = find_openmeeg_include()
    #openmeeg_lib = (openmeeg_root / '..' / 'lib').resolve()
    #assert openmeeg_include.is_dir()
    #assert (openmeeg_include / 'OpenMEEG' / 'vect3.h').is_file()
    #numpy_include = np.get_include()
    #swig_openmeeg = Extension(
    #    "openmeeg._openmeeg",
    #    ["openmeeg/openmeeg.i", "openmeeg/openmeeg.cpp"],
    #    include_dirs=[openmeeg_include],
    #    swig_opts=['-c++', '-v', '-Werror', f'-I{openmeeg_include}'],
    #    libraries=['OpenMEEG'],
    #    library_dirs=[openmeeg_lib],
    #    extra_compile_args=['-v', '-std=c++17'],
    #)

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
          #ext_modules=[swig_openmeeg],
    )
