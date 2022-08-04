|GitHub Actions|_ |lgtm|_ |CodeCov|_ |PyPIVersion|_ |condaVersion|_

.. |GitHub Actions| image:: https://github.com/openmeeg/openmeeg/actions/workflows/build_and_test.yml/badge.svg
.. _Github Actions: https://github.com/openmeeg/openmeeg/actions/workflows/build_and_test.yml

.. |lgtm| image:: https://img.shields.io/lgtm/grade/cpp/g/openmeeg/openmeeg.svg?logo=lgtm&logoWidth=18
.. _lgtm: https://lgtm.com/projects/g/openmeeg/openmeeg/context:cpp

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

Unix (Linux & Mac OS X)
^^^^^^^^^^^^^^^^^^^^^^^

On Debian/Ubuntu you will need to install the dependencies with::

    $ sudo apt-get install gcc g++ make cmake libopenblas-dev liblapacke-dev libmatio-dev libhdf5-dev

*optionally*::

    $ sudo apt-get install python3-numpy swig libvtk6-dev doxygen graphviz libcgal-dev

On Fedora and Centos::

    $ sudo yum install gcc make cmake openblas-devel hdf5-devel matio-devel

*optionally*::

    $ sudo yum install python3-numpy swig vtk-devel doxygen cgal-devel

To install OpenMEEG from source in a terminal::

    $ git clone https://github.com/openmeeg/openmeeg.git

then::

    $ cd openmeeg
    $ cmake -D build -DCMAKE_BUILD_TYPE=Release -DUSE_PROGRESSBAR=ON -DBLA_VENDOR=OpenBLAS .
    $ cmake --build build --config=Release

**Note for Python users**:

- To use Python bindings you will need a recent version of CMake >= 3.16.2
- and a recent version of Swig >= 4.0

You will need to define more CMake variables if you want the support for:

`-DENABLE_PYTHON=ON`` (Python >= 3.7 is required)
    Enable Python wrapping.
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

You can run the full test suite with::

    $ pushd build
    $ ctest -C Release || ctest -C Release --rerun-failed --output-on-failure
    $ popd

If no test is failing you can install with (and optionally with ``--install-prefix=...`` to install somewhere other than the default)::

    $ cmake --build build --target install

You should now be able to run the *om_assemble* command and see something like this::

    $ om_assemble
    om_assemble version 2.4.7 compiled at Jul 26 2022 18:17:12

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

Windows
^^^^^^^

You will need to install MSVC 15 (2017) or later and `CMake <http://www.cmake.org>`_,
which can be installed via ``pip``.
Then download the source from github, and follow the steps that we use to
build OpenMEEG on GitHub Actions: `.github/workflows/build_and_test.yml <https://github.com/openmeeg/openmeeg/blob/main/.github/workflows/build_and_test.yml>`_

Supported Blas/Lapack Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We support `OpenBLAS <http://www.openblas.net/>`_ and
`Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_
on Linux, macOS, and Windows; and `ATLAS <http://math-atlas.sourceforge.net/>`_
and reference Netlib `BLAS <https://netlib.org/blas/>`_ /
`LAPACK <https://netlib.org/lapack/>`_ on Linux.

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
