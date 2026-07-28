|GitHub Actions|_ |CodeQL|_ |CodeCov|_ |PyPIVersion|_ |condaVersion|_

.. |GitHub Actions| image:: https://github.com/openmeeg/openmeeg/actions/workflows/build_and_test.yml/badge.svg
.. _Github Actions: https://github.com/openmeeg/openmeeg/actions/workflows/build_and_test.yml

.. |CodeQL| image:: https://github.com/openmeeg/openmeeg/workflows/CodeQL/badge.svg
.. _CodeQL: https://github.com/openmeeg/openmeeg/actions/workflows/codeql-analysis.yml

.. |CodeCov| image:: https://codecov.io/gh/openmeeg/openmeeg/branch/main/graph/badge.svg
.. _CodeCov: https://codecov.io/gh/openmeeg/openmeeg

.. |PyPIVersion| image:: https://badge.fury.io/py/openmeeg.svg
.. _PyPIVersion: https://badge.fury.io/py/openmeeg

.. |condaVersion| image:: https://anaconda.org/conda-forge/openmeeg/badges/version.svg
.. _condaVersion: https://anaconda.org/conda-forge/openmeeg

OpenMEEG: forward problems solver in the field of EEG and MEG
=============================================================

.. highlight:: console

The OpenMEEG software is a C++ package for solving the forward
problems of electroencephalography (EEG) and magnetoencephalography (MEG).

OpenMEEG is distributed under the French opensource license CeCILL-B. It is
intended to give users the freedom to modify and redistribute the software.
It is therefore compatible with popular opensource licenses such as the GPL
and BSD licenses. The CeCILL-B license imposes to anybody distributing a
software incorporating OpenMEEG the obligation to give credits (by citing the
appropriate publications), in order for all contributions to be properly
identified and acknowledged.

Cite this software
------------------

The references to be acknowledged are ::

    Gramfort A, Papadopoulo T, Olivi E, Clerc M. OpenMEEG: opensource software for quasistatic
    bioelectromagnetics. Biomedical engineering online (2010) vol. 9 (1) pp. 45

    Kybic J, Clerc M, Abboud T, Faugeras O, Keriven R, Papadopoulo T. Generalized head models for MEG/EEG: boundary element method
    beyond nested volumes. Phys. Med. Biol. (2006) vol. 51 pp. 1333-1346

.. image:: https://raw.githubusercontent.com/openmeeg/openmeeg.github.io/source/_static/inria.png

Install precompiled library and Python bindings
-----------------------------------------------

To install OpenMEEG (along with the binary applications) via `anaconda <https://www.anaconda.com/download/>`_ you can just do::

    $ conda install -c conda-forge openmeeg

Python wrappers can also be installed via `pip`::

    $ pip install openmeeg

On Fedora::

    $ dnf install openmeeg openmeeg-devel python2-openmeeg

On RHEL/CentOS 7, enable `EPEL repositories <https://fedoraproject.org/wiki/EPEL>`_ and install::

    $ yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
    $ yum install openmeeg openmeeg-devel python2-openmeeg

Additional repositories recommended on RHEL 7::

    $ subscription-manager repos --enable "rhel-*-optional-rpms" --enable "rhel-*-extras-rpms"

Build OpenMEEG from source
--------------------------

On any operating system, you should get the latest OpenMEEG source the usual way::

    $ git clone https://github.com/openmeeg/openmeeg.git
    $ cd openmeeg

Then you need to get dependencies installed and configured for your operating system.

Building on Linux
^^^^^^^^^^^^^^^^^

On Debian/Ubuntu you will need to install the dependencies with (Fedora flavors can use a similar command)::

    $ sudo apt install gcc g++ make cmake libopenblas-dev liblapacke-dev libmatio-dev libhdf5-dev

*optionally*::

    $ sudo apt install python3-numpy swig libvtk6-dev doxygen graphviz libcgal-dev

then::

    $ ./build_tools/cmake_configure.sh
    $ cmake --build build --config release

The ``cmake_configure.sh`` script should automatically set the build to configure
Python using SWIG.

Building on macOS
^^^^^^^^^^^^^^^^^
For local debugging, it's easiest to use ``brew`` to install dependencies::

    $ brew install hdf5 libmatio libomp swig openblas

Then follow brew's suggestion to add to your paths (probably in ``.bash_profile`` or some similar place) with something like the following::

    $ export PATH="$HOMEBREW_PREFIX/opt/llvm/bin:$PATH
    $ export LDFLAGS="-L$HOMEBREW_PREFIX/opt/llvm/lib -L$HOMEBREW_PREFIX/opt/openblas/lib"
    $ export CPPFLAGS="-I$HOMEBREW_PREFIX/opt/llvm/include -I$HOMEBREW_PREFIX/opt/openblas/include"

Then you should be able to build as usual::

    $ ./build_tools/cmake_configure.sh
    $ cmake --build build --config release

Building on Windows
^^^^^^^^^^^^^^^^^^^
One configuration that makes Windows development easier is getting a usable
Bash shell under Windows properly configured to compile using Visual Studio.
The steps are roughly:

1. Install some variant of `Visual Studio <https://visualstudio.microsoft.com/downloads/>`__ (e.g., `2019 <https://visualstudio.microsoft.com/vs/older-downloads/>`__), the community variants are free and should work.
2. Install the `Git for Windows SDK <https://github.com/git-for-windows/build-extra/releases>`_.
3. Launch a ``x64 Native Tools Command Prompt for VS 2019`` (i.e., a variant of ``cmd``),
   which can be done from the Start menu.
4. Run ``C:\git-sdk-64\usr\bin\bash -l`` from within that prompt.

.. note:: If you do not have access to Windows but need to debug it, consider
          using the `Windows VM dev images <https://developer.microsoft.com/en-us/windows/downloads/virtual-machines/>`__.

For dependencies on Windows, we make use of ``vcpkg``. The default generator
is ``"Visual Studio 17 2022"``, if you would like to use 2019 then set::

    $ export CMAKE_GENERATOR="Visual Studio 16 2019"

Then you can use our convenience script for setting up ``vcpkg``. From the ``openmeeg``
root, run::

    $ source ./build_tools/setup_vcpkg_compilation.sh

Then you need MKL or OpenBLAS. The easiest way to get this is to use ``brew`` on macOS,
``apt`` or similar on Linux, and download precompiled OpenBLAS binaries on Windows
from using ``tools/download_openblas_windows.sh``::

    $ ./build_tools/download_openblas_windows.sh
    $ export CMAKE_PREFIX_PATH=$(cygpath -m $PWD/openblas/64)
    $ export CMAKE_CXX_FLAGS="-I$(cygpath -m $PWD/openblas/64/include)"
    $ export PATH=$PATH:$PWD/openblas/64/lib

Then you also need the path to the compiled libraries for tests to work::

    $ export PATH=$PATH:$PWD/build/OpenMEEG/Release:$PWD/build/OpenMEEGMaths/Release

.. note:: Consider adding ``export`` statements to your ``~.bashrc`` to
          facilitate future debugging, but be sure to translate the ``$PWD``
          to the actual Unix-formatted path on your system, e.g.::

              export CMAKE_GENERATOR="Visual Studio 16 2019"
              export CMAKE_PREFIX_PATH=C:/Users/whoever/python/openmeeg/openblas/64
              export CMAKE_CXX_FLAGS="-IC:/Users/whoever/python/openmeeg/openblas/64/include"
              export PATH=$PATH:/c/Users/whoever/python/openmeeg/openblas/64/lib:/c/Users/whoever/python/openmeeg/build/OpenMEEG/Release:/c/Users/whoever/python/openmeeg/build/OpenMEEGMaths/Release

Then you can build as usual::

    $ ./build_tools/cmake_configure.sh
    $ cmake --build build --config release

The configure step will take a few minutes because this is the stage during
which ``vcpkg`` builds dependencies (and HDF5 in particular takes some time).
But once it has completed, any subsequent ``./build_tools/cmake_configure.sh``
calls should be much faster because the completed dependency builds are stored
in the ``vcpkg`` directory for future use.

If you ever have problems with DLL linkage, consider using::

    $ ./build_tools/install_dependency_walker.sh
    $ ./Dependencies/DependenciesGui.exe

to examine issues with ``OpenMEEGMaths.dll`` for example.

Testing
^^^^^^^
Once you have a complete build in ``build``, you can test with::

    $ cd build
    $ ctest -C release || ctest -C release --rerun-failed --output-on-failure

Optional build variables
^^^^^^^^^^^^^^^^^^^^^^^^
You will need to define more CMake variables if you want the support for:

`-DENABLE_PYTHON=ON`` (Python >= 3.10 is required)
    Enable Python wrapping (automatically enabled by cmake_configure.sh)
`-DUSE_VTK=ON`
    VTK file format support.
`-DUSE_CGAL=ON`
    CGAL meshing tools.
`-DBUILD_DOCUMENTATION=ON`
    Reference documentation. Make sure to have `doxygen` with `dot` support.
`-DENABLE_WERROR=ON`
    Treat compilation warnings as errors
`-DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache`
    To speed up builds. `ccache` must be installed.

Installation
^^^^^^^^^^^^
In usual cmake fashion, you can install with (and optionally with ``--install-prefix=...`` to install somewhere other than the default)::

    $ cmake --build build --target install

You should now be able to run the *om_assemble* command and see something like this::

    $ om_assemble
    om_assemble version 2.5.5 compiled at Aug 26 2022 18:17:12

    om_assemble [-option] [filepaths...]

    option :
       -HeadMat, -HM, -hm :
           Compute Head Matrix for Symmetric BEM (left-hand side of linear system).
           ...

In some Linux distributions (AMD64/X86_64) you may see some errors like this::

    Error while loading shared libraries: libOpenMEEG.so.1: cannot open shared object file: No such file or directory

You need to ensure that the ``install`` target libraries (given the prefix that
was used) is in your library search path, e.g., by settincg ``LD_LIBRARY_PATH``
or editing ``/etc/ld.so.conf`` and using ``sudo ldconfig``.

You can now give a try to OpenMEEG on the `sample dataset <https://github.com/openmeeg/openmeeg_sample_data/archive/master.zip>`_.

Supported Blas/Lapack Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We support `OpenBLAS <http://www.openblas.net/>`_ and
`Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ on Linux, macOS, and Windows.

Using OpenMEEG
--------------

Have a look into the `tutorial <https://openmeeg.github.io/tutorial.html>`_
for more info and for defining your geometry.

CeCILL-B full license
---------------------

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's authors, the holders of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
