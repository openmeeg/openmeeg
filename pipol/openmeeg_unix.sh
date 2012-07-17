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

sh ~/install_packages.sh

BRANCH=master

git clone --recursive git://github.com/openmeeg/openmeeg.git

perl ./openmeeg/pipol/cmake.pl
cd openmeeg

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
    limit stacksize unlimited # For integrated lapack testing
    ctest -D ExperimentalTest
    ctest -D ExperimentalSubmit
    make package
    make clean
    mv OpenMEEG*tar.gz ~/.pipol/packages/openmeeg-${BRANCH}
    mv OpenMEEG*dmg* ~/.pipol/packages/openmeeg-${BRANCH}
    mv OpenMEEG*rpm* ~/.pipol/packages/openmeeg-${BRANCH}
    mv OpenMEEG*deb* ~/.pipol/packages/openmeeg-${BRANCH}
}

OMP=OFF
if [ x$SYSTEM = xLinux ] ; then
    if [ x$ARCH = xx86_64 ]; then
        OMP=ON
    fi
fi

# Blas Lapack (Atlas)
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DENABLE_PYTHON=ON -DUSE_PROGRESSBAR=ON -DUSE_OMP=${OMP} .
build

# MKL (static mode)
rm -f CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DUSE_ATLAS=OFF -DUSE_MKL=ON -DBUILD_SHARED_LIBS=OFF -DUSE_PROGRESSBAR=ON -DUSE_OMP=${OMP} .
build

# Pure blas/Lapack (No atlas)
rm -f CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DENABLE_PYTHON=ON -DUSE_ATLAS=OFF -DUSE_MKL=OFF -DUSE_PROGRESSBAR=ON -DUSE_OMP=${OMP} .
build
