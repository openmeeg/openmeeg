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

.. image:: http://biomaj.genouest.org/wp-content/uploads/2011/07/logo-inria_im.png

Install precompiled binaries
----------------------------

`Download precompiled binaries <https://gforge.inria.fr/frs/?group_id=435>`_ of the latest stable version.

Or on Mac OS install with `Homebrew <http://mxcl.github.com/homebrew/>`_::

    $ brew install openmeeg

On Ubuntu/Debian GNU Linux you can use the http://neuro.debian.net package repository.

Build OpenMEEG
--------------

Unix (Linux & Mac OS X)
^^^^^^^^^^^^^^^^^^^^^^^

On Ubuntu/Debian you will need to install the dependencies with::

    sudo apt-get install git-core gcc g++ make cmake libatlas-base-dev python-numpy python-dev swig

On Fedora and Centos::

    sudo yum install git-core gcc make cmake wget perl atlas-devel blas-devel numpy python-devel

On Mac OS X, you'll need `CMake <http://www.cmake.org>`_ (install it with `Homebrew <http://mxcl.github.com/homebrew/>`_ or `Fink <http://www.finkproject.org/>`_ or `Macports <http://www.macports.org/>`_ or by direct download).

Then from a terminal::

    $ git clone --recursive git://github.com/openmeeg/openmeeg.git

or if it does not work try::

    $ git clone --recursive https://github.com/openmeeg/openmeeg.git

then::

    $ cd openmeeg
    $ mkdir build
    $ cd build
    $ cmake -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DUSE_PROGRESSBAR=ON ..
    $ make

If you want the support for:
-Non-nested geometries, you will need to add "-DUSE_VTK=ON" to the cmake line above.
-Python, you will need to add "-DENABLE_PYTHON=ON".
-parallel computation with OpenMP, add "-DUSE_OMP=ON".

Then you can run the test suite with::

    $ make test

If no test is failing you can install with::

    $ make install

You should now be able to run the *om_assemble* command and see something like this::

    $ om_assemble
    om_assemble version 2.2.dev (802) compiled at Sep 20 2011 11:50:08

    Not enough arguments
    Please try "om_assemble -h" or "om_assemble --help "

You can now give a try to OpenMEEG on the `sample dataset <https://gforge.inria.fr/frs/download.php/29059/openmeeg_sample_dataset.zip>`_.

Windows
^^^^^^^

You will need to install visual studio, `CMake <http://www.cmake.org>`_ and the
`Intel MKL library <http://software.intel.com/en-us/intel-mkl/>`_.
Then download the source from github, load the CMake.exe GUI, set the proper option
and generate the visual studio project. You can then open it and build the project.

Using OpenMEEG
--------------

Have a look into data/README.rst for defining your geometry.

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
