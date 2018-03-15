#! /bin/sh

basedir=`dirname $0`

if [ x$1 != "x--incremental" ]; then
    rm -rf build
    mkdir build
fi

cd build
cmake -DENABLE_PYTHON=ON -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DMATLAB_TESTING=OFF ..
ctest -D ExperimentalConfigure
ctest -D ExperimentalBuild
ctest -D ExperimentalTest --no-compress-output || true

if [ -f Testing/TAG ]; then
    xsltproc ../${basedir}/CTest2JUnit.xsl Testing/`head -n 1 < Testing/TAG`/Test.xml > JUnitTestResults.xml
fi

# cdash backward compatibility.

ctest -D ExperimentalSubmit
