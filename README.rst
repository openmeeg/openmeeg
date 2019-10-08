|Travis|_ |AppVeyor|_ |lgtm|_ |CodeCov|_ |condaVersion|_ |gitter|_

.. |Travis| image:: https://api.travis-ci.org/openmeeg/openmeeg.svg?branch=master
.. _Travis: https://travis-ci.org/openmeeg/openmeeg/branches

.. |AppVeyor| image:: https://ci.appveyor.com/api/projects/status/11um4d4c8nn4itju/branch/master?svg=true
.. _AppVeyor: https://ci.appveyor.com/project/openmeegci/openmeeg/history

.. |CodeCov| image:: https://codecov.io/gh/openmeeg/openmeeg/branch/master/graph/badge.svg
.. _CodeCov: https://codecov.io/gh/openmeeg/openmeeg

.. |condaVersion| image:: https://anaconda.org/conda-forge/openmeeg/badges/version.svg
.. _condaVersion: https://anaconda.org/conda-forge/openmeeg

.. |gitter| image:: https://badges.gitter.im/openmeeg/openmeeg.svg
.. _gitter: https://gitter.im/openmeeg/openmeeg

.. |lgtm| image:: https://img.shields.io/lgtm/grade/cpp/g/openmeeg/openmeeg.svg?logo=lgtm&logoWidth=18
.. _lgtm: https://lgtm.com/projects/g/openmeeg/openmeeg/context:cpp

OpenMEEG: forward problems solver in the field of EEG and MEG
=============================================================

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

Install precompiled binaries
----------------------------

Binaries for Linux/Mac/Windows are available at `Download precompiled binaries <http://openmeeg.gforge.inria.fr/download/?C=M;O=D>`_.

To install OpenMEEG via `anaconda <https://www.anaconda.com/download/>`_ you can just do::

    $ conda install -c conda-forge openmeeg


On Ubuntu/Debian GNU Linux you may be able use the http://neuro.debian.net package repository.

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

    sudo apt-get install gcc g++ make cmake libopenblas-dev liblapacke-dev libhdf5-serial-dev libmatio-dev

*optionally*::

    sudo apt-get install python-dev python-numpy swig libvtk6-dev doxygen graphviz libcgal-dev

On Fedora and Centos::

    sudo yum install gcc make cmake openblas-devel hdf5-devel matio-devel

*optionally*::

    sudo yum install python-devel python2-numpy swig vtk-devel doxygen cgal-devel

On Mac OS X, we recommend you install it with `Homebrew <http://brew.sh/>`_ using::

    $ brew install https://raw.githubusercontent.com/openmeeg/openmeeg/master/openmeeg.rb

To install with Homebrew the current development version use::

    $ brew install https://raw.githubusercontent.com/openmeeg/openmeeg/master/openmeeg.rb --devel

To install OpenMEEG from source in a terminal::

    $ git clone https://github.com/openmeeg/openmeeg.git

then::

    $ cd openmeeg
    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PROGRESSBAR=ON -DBLA_VENDOR=OpenBLAS ..
    $ make


You will need to define more CMake variables if you want the support for:

- Python wrapping, add "-DENABLE_PYTHON=ON" (default version is 3, add "-DPYTHON_VERSION=2" to use python2)

- VTK file format, add "-DUSE_VTK=ON".

- CGAL meshing tools, add "-DUSE_CGAL=ON".

- Reference documentation, add "-DBUILD_DOCUMENTATION=ON". Make sure to have `doxygen` with `dot` support.

Then you can run the full test suite with::

    $ make

or if you just want to run the tests for OpenMEEG::

    $ make test

If no test is failing you can install with::

    $ make install

You should now be able to run the *om_assemble* command and see something like this::

    $ om_assemble
    om_assemble version 2.4.0 compiled at Mar 21 2018 18:17:12

    om_assemble [-option] [filepaths...]

    option :
       -HeadMat, -HM, -hm :
           Compute Head Matrix for Symmetric BEM (left-hand side of linear system).
           ...

In some Linux distributions (AMD64/X86_64) you may see some errors like this::

    Error while loading shared libraries: libOpenMEEG.so.1: cannot open shared object file: No such file or directory

OpenMEEG puts its libraries in "/usr/local/lib64", which is not included
in your loader's search path. If so, run this command as root::

    # echo '/usr/local/lib64/' >> /etc/ld.so.conf && ldconfig

Now you can try to run the *om_assemble* again.

You can now give a try to OpenMEEG on the `sample dataset <https://github.com/openmeeg/openmeeg_sample_data/archive/master.zip>`_.

Windows
^^^^^^^

You will need to install visual studio, `CMake <http://www.cmake.org>`_.
Then download the source from github, load the CMake.exe GUI, set the proper option
and generate the visual studio project. You can then open it and build the project.
Note that on Windows we currently recommend to use Intel MKL library.
See how we build OpenMEEG on AppVeyor: `.appveyor.yml <https://github.com/openmeeg/openmeeg/blob/master/.appveyor.yml>`_

Supported Blas-Lapack Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- on Linux: `Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ , `OpenBLAS <http://www.openblas.net/>`_ (and possibly `Atlas <http://math-atlas.sourceforge.net>`_)

- on Mac OS X: `Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ , `OpenBLAS <http://www.openblas.net/>`_, `vecLib <https://developer.apple.com/reference/accelerate/veclib>`_

- on Windows: `Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ , `OpenBLAS <http://www.openblas.net/>`_

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
