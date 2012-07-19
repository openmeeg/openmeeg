#!/bin/bash

SYSTEM=`uname`
ARCH=`uname -m`
HOME=`pwd`

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

#   Attempt to determine package type.

PACKAGE_TYPE=
PACKAGING_OPTION=
if [ x$SYSTEM = xLinux ] ; then
    if [ -e /usr/bin/yum ] ; then
        PACKAGE_TYPE=rpm
        PACKAGING_OPTION="-DBUILD_RPM=ON"
    else
        if [ -e /usr/bin/apt-get ] ; then
            PACKAGE_TYPE=deb
            #PACKAGING_OPTION="-DBUILD_DEB=ON"
        else
            PACKAGE_TYPE=dmg
        fi
    fi
fi

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

# Needed for lapack tests.

ulimit -s unlimited

function build {
    ctest -D ExperimentalConfigure
    ctest -D ExperimentalBuild
    ctest -D ExperimentalTest
    ctest -D ExperimentalSubmit
    make package
    mv OpenMEEG*.tar.* ~/.pipol/packages/openmeeg-${BRANCH}
    echo "OpenMEEG*.tar.*"
    ls -l ~/.pipol/packages/openmeeg-${BRANCH}
    if [ x$PACKAGE_TYPE != x ] ; then
        make OpenMEEG_${PACKAGE_TYPE}
        mv OpenMEEG*${PACKAGE_TYPE}* ~/.pipol/packages/openmeeg-${BRANCH}
    fi
    mkdir -p ~/.pipol/$PIPOL_HOST/
    file=`mktemp -d --tmpdir=${HOME}/.pipol/$PIPOL_HOST/ tmpXXXX`
    cp contrib/matio/test/MATIO* $file
    make clean
}

OMP=OFF
if [ x$SYSTEM = xLinux ] ; then
    if [ x$ARCH = xx86_64 ]; then
        OMP=ON
    fi
fi

# Blas Lapack (Atlas)
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON ${PACKAGING_OPTION} -DENABLE_PYTHON=ON -DUSE_PROGRESSBAR=ON -DUSE_OMP=${OMP} .
build

# MKL (static mode)
rm -f CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON ${PACKAGING_OPTION} -DUSE_ATLAS=OFF -DUSE_MKL=ON -DBUILD_SHARED_LIBS=OFF -DUSE_PROGRESSBAR=ON -DUSE_OMP=${OMP} .
build

# Pure blas/Lapack (No atlas)
rm -f CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON ${PACKAGING_OPTION} -DENABLE_PYTHON=ON -DUSE_ATLAS=OFF -DUSE_MKL=OFF -DUSE_PROGRESSBAR=ON -DUSE_OMP=${OMP} .
build
