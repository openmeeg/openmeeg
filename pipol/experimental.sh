#!/bin/bash

SYSTEM=`uname`

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

cmake -DBUILD_TESTING=True -DENABLE_PACKAGING=True -DPYTHON_WRAP=True .
# cmake -DBUILD_TESTING=True -DENABLE_PACKAGING=True -DBUILD_SHARED=False -DCMAKE_CXX_FLAGS="-fPIC"  -DCMAKE_C_FLAGS="-fPIC" .
ctest -D ExperimentalStart
ctest -D ExperimentalConfigure
ctest -D ExperimentalBuild
ctest -D ExperimentalTest
ctest -D ExperimentalSubmit
make clean
