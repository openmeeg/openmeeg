#!/usr/bin/env python

# Copyright (C) 2011-2020 Alexandre Gramfort
# <alexandre.gramfort@inria.fr>
# This assumes/requires that SWIG has already been run so the .so exists!

from pathlib import Path
import os
import sys

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


if __name__ == "__main__":
    import numpy as np
    manifest = (root / 'MANIFEST')
    if manifest.is_file():
        os.remove(manifest)

    with open('README.rst', 'r') as fid:
        long_description = fid.read()

    # SWIG
    cmdclass = {}
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
            library_dirs.append(str(openmeeg_lib))
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
    else:  # built with -DENABLE_PYTHON=ON
        # TODO: This breaks macOS for some reason!
        if sys.platform != 'darwin':
            cmdclass['bdist_wheel'] = bdist_wheel

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
          cmdclass=cmdclass,
          ext_modules=ext_modules,
          )
