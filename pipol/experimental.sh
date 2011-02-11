#!/bin/bash

SYSTEM=`uname`
ARCH=`uname -m`

if [ -e ./pipol ] ; then
	rm -rf ./pipol/$PIPOL_HOST
	mkdir ./pipol/$PIPOL_HOST
else
	mkdir ./pipol
	rm -rf ./pipol/$PIPOL_HOST
	mkdir ./pipol/$PIPOL_HOST
fi
cd ./pipol/$PIPOL_HOST

svn checkout svn://scm.gforge.inria.fr/svn/openmeeg/trunk openmeeg-trunk --quiet

sh ./openmeeg-trunk/pipol/install_packages.sh
perl ./openmeeg-trunk/pipol/cmake.pl

cd openmeeg-trunk

# Handle MKL
if [ x$SYSTEM = xDarwin ] ; then
    export MKLDIR=/netshare/i386_mac/icc/11.0.074/Frameworks/mkl/ # mac 32
fi
if [ x$SYSTEM = xLinux ] ; then
    if [ x$ARCH = xx86_64 ]; then
        export MKLDIR=/netshare/amd64/icc/11.0.074/mkl/ # linux 64
    else
        export MKLDIR=/netshare/i386/icc/11.0.074/mkl/ # linux 32
    fi
fi

function build {
    ctest -D ExperimentalConfigure
    ctest -D ExperimentalBuild
    ctest -D ExperimentalTest
    ctest -D ExperimentalSubmit
    make package
    make clean
    mv OpenMEEG*tar.gz ~/.pipol/packages
    mv OpenMEEG*dmg* ~/.pipol/packages
    mv OpenMEEG*rpm* ~/.pipol/packages
    mv OpenMEEG*deb* ~/.pipol/packages
}

# Blas Lapack
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DENABLE_PYTHON=ON -DUSE_PROGRESSBAR=ON .
build

# MKL (static mode)
rm -f CMakeCache.txt
# cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DENABLE_PYTHON=ON -DUSE_ATLAS=OFF -DUSE_MKL=ON -DBUILD_SHARED=OFF -DUSE_PROGRESSBAR=ON .
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DUSE_ATLAS=OFF -DUSE_MKL=ON -DBUILD_SHARED=OFF -DUSE_PROGRESSBAR=ON .
build

