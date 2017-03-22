|Travis|_ |AppVeyor|_ |CodeCov|_

.. |Travis| image:: https://api.travis-ci.org/openmeeg/openmeeg.svg?branch=master
.. _Travis: https://travis-ci.org/openmeeg/openmeeg

.. |AppVeyor| image:: https://ci.appveyor.com/api/projects/status/github/openmeeg/openmeeg?branch=master&svg=true
.. _AppVeyor: https://ci.appveyor.com/project/agramfort/openmeeg/history

.. |CodeCov| image:: https://codecov.io/gh/openmeeg/openmeeg/branch/master/graph/badge.svg
.. _CodeCov: https://codecov.io/gh/openmeeg/openmeeg


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

    Gramfort et al. OpenMEEG: opensource software for quasistatic
    bioelectromagnetics. Biomedical engineering online (2010) vol. 9 (1) pp. 45

    Kybic et al. Generalized head models for MEG/EEG: boundary element method
    beyond nested volumes. Phys. Med. Biol. (2006) vol. 51 pp. 1333-1346

.. image:: https://raw.githubusercontent.com/openmeeg/openmeeg.github.io/source/_static/inria.png

Install precompiled binaries
----------------------------

`Download precompiled binaries <http://openmeeg.gforge.inria.fr/download/>`_

On Ubuntu/Debian GNU Linux you can use the http://neuro.debian.net package repository.

Build OpenMEEG
--------------

Unix (Linux & Mac OS X)
^^^^^^^^^^^^^^^^^^^^^^^

On Debian/Ubuntu you will need to install the dependencies with::

    sudo apt-get install git-core gcc g++ make cmake libopenblas-dev liblapacke-dev

*optionally*::

    sudo apt-get install libhdf5-serial-dev libmatio-dev python-dev python-numpy swig libvtk6-dev doxygen libcgal-dev

On Fedora and Centos::

    sudo yum install git-core gcc make cmake wget perl openblas-devel

*optionally*::

    sudo yum install hdf5-devel matio-devel python-devel python2-numpy swig vtk-devel doxygen cgal-devel

On Mac OS X, you'll need `CMake <http://www.cmake.org>`_ (install it with `Homebrew <http://brew.sh/>`_ or `Fink <http://www.finkproject.org/>`_ or `Macports <http://www.macports.org/>`_ or by direct download).

    *e.g* with Homebrew::

    $ brew tap homebrew/science
    $ brew tap homebrew/python
    $ brew update && brew upgrade

    *optionally*::

    $ brew install hdf5 && brew install libmatio --with-hdf5
    $ brew install openblas && brew link openblas --force
    $ brew install vtk
    $ brew install --with-qt5 cgal
    $ brew install Doxygen

Then from a terminal::

    $ git clone git://github.com/openmeeg/openmeeg.git

or if it does not work try::

    $ git clone https://github.com/openmeeg/openmeeg.git

then::

    $ cd openmeeg
    $ mkdir build
    $ cd build
    $ cmake -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DUSE_PROGRESSBAR=ON -DBUILD_DOCUMENTATION=OFF -DENABLE_PYTHON=OFF ..
    $ make

OpenMEEG will download and compile the **zlib**, (optional **vtk**), **hdf5**, and **matio** by default. In case your system already provides these libraries (see optional packages), you should specify the following variables to the cmake command line: "-DUSE_SYSTEM_zlib=ON -DUSE_SYSTEM_hdf5=ON -DUSE_SYSTEM_matio=ON".

You will need to define more CMake variables if you want the support for:

- Parallel computation with OpenMP, add "-DUSE_OMP=ON". Then you can set the number of threads you want to run in parallel by setting the OMP_NUM_THREADS environment variable. On a Unix system if you want to run 4 threads in parallel do (before running OpenMEEG)::

    $ export OMP_NUM_THREADS=4

- Python wrapping, add "-DENABLE_PYTHON=ON".

- VTK file format, add "-DUSE_VTK=ON".

- CGAL meshing tools, add "-DUSE_CGAL=ON".

Then you can run the full test suite with::

    $ make check

or if you just want to run the tests for OpenMEEG::

    $ make test-OpenMEEG

If no test is failing you can install with::

    $ make install

You should now be able to run the *om_assemble* command and see something like this::

    $ om_assemble
    om_assemble version 2.3.dev compiled at Oct 30 2016 18:43:25

    Not enough arguments
    Please try "om_assemble -h" or "om_assemble --help "

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

Supported Blas-Lapack Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- on Linux: `Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ , `OpenBLAS <http://www.openblas.net/>`_, `Atlas <http://math-atlas.sourceforge.net>`_, Lapack

- on Mac OS X: `Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ , `OpenBLAS <http://www.openblas.net/>`_, `vecLib <https://developer.apple.com/reference/accelerate/veclib>`_

- on Windows: `Intel MKL <http://software.intel.com/en-us/intel-mkl/>`_ , `Clapack <https://github.com/openmeeg/clapack>`_

Using OpenMEEG
--------------

Have a look into data/README.rst for defining your geometry.
and/or
specify to cmake "-DBUILD_TUTORIALS=ON", and read the tutorials in pdf.

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
