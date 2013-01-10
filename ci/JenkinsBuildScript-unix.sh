#! /bin/sh

basedir=`dirname $0`

if [ x$1 != "x--incremental" ]; then
    rm -rf build
    mkdir build
fi

cd build
cmake -DUSE_ATLAS=ON -DUSE_MKL=OFF -DENABLE_PYTHON=ON -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DMATLAB_TESTING=OFF ..
ctest -D ExperimentalConfigure
ctest -D ExperimentalBuild
ctest -D ExperimentalTest
ctest -D ExperimentalSubmit

xsltproc ${basedir}/ci/ctest2junix.xsl Testing/`head -n 1 < Testing/TAG`/Test.xml > JUnitTestResults.xml
