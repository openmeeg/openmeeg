#!/bin/bash

#PIPOL esn amd64_kvm-linux-debian-lenny  none 2:00 --silent --root
#PIPOL esn amd64-linux-fedora-core11.dd.gz none 2:00 --silent --root
#PIPOL esn i386_kvm-linux-debian-lenny  none 2:00 --silent --root
#PIPOL esn i386-linux-fedora-core11.dd.gz none 2:00 --silent --root
#PIPOL esn x86_64_mac-mac-osx-server-snow-leopard none 0:40 --silent --root
#PIPOL esn i386_mac-mac-osx-server-leopard none 00:40 --silent --root

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
    make NightlyUpdate
    make NightlyConfigure
    make NightlyBuild
    make NightlyTest
    if [ x$SYSTEM = xLinux ] ; then
        if [ x$ARCH = xx86_64 ]; then
            make NightlyCoverage
            make NightlyMemCheck
        fi
    fi
    make NightlySubmit
    make clean
}


if [ x$SYSTEM = xLinux ] ; then
    if [ x$ARCH = xx86_64 ]; then
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DENABLE_COVERAGE=ON .
    else
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DENABLE_PYTHON=ON .
    fi
fi

if [ x$SYSTEM = xDarwin ] ; then
    cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON .
fi

build

## MKL (static mode)
#rm -f CMakeCache.txt
#cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DUSE_ATLAS=OFF -DUSE_MKL=ON -DBUILD_SHARED=OFF -DUSE_PROGRESSBAR=ON .
#build

