#!/bin/bash

#PIPOL esn i386_kvm-windows-xp-pro-sp3 none 02:00 --silent --root
#PIPOL esn amd64_kvm-windows-7 none 02:00 --silent --root

function build_and_test() {
    # 2) Creation of temporary repertory where to build
    export TEMP=/cygdrive/c/Temp

    rm -rf ${TEMP}
    mkdir ${TEMP}

    # 3) Configuration of the variables for matlab
    if expr $PIPOL_IMAGE : ".*amd64.*"; then
        SYSTEMOS=" Win64"
        SYSTEMCOMPILE="x64"
        export PROCESSOR_ARCHITECTURE=AMD64
        export MKL_DIR=Y:/amd64/icc-windows-2008/Compiler/11.1/046/mkl
        export ICC_LIB_DIR=Y:amd64/icc-windows-2008/Compiler/11.1/046/lib/intel64/
        export MKL_LIB_DIR=Y:amd64/icc-windows-2008/Compiler/11.1/046/mkl/em64t/lib
    else
        SYSTEMCONF=""
        SYSTEMCOMPILE="win32"
        export PROCESSOR_ARCHITECTURE=x86
        export MKL_DIR=Y:i386/icc-windows-xp/Compiler/11.1/046/mkl/
        export ICC_LIB_DIR=Y:i386/icc-windows-xp/Compiler/11.1/046/lib/ia32/
        export MKL_LIB_DIR=Y:i386/icc-windows-xp/Compiler/11.1/046/mkl/ia32/lib/
    fi

    # 4) Detection of cmake version

    wget http://www.cmake.org/files/v2.8/cmake-2.8.8-win32-x86.zip
    unzip cmake-2.8.8-win32-x86.zip
    CMAKE = `pwd`/cmake-2.8.8-win32-x86/bin/cmake.exe
    CPACK = `pwd`/cmake-2.8.8-win32-x86/bin/cpack.exe
    CTEST = `pwd`/cmake-2.8.8-win32-x86/bin/ctest.exe

    # 5) Cleaning the projet
    cd ${TEMP}

    BRANCH=master
    git clone --recursive git://github.com/openmeeg/openmeeg.git
    cd openmeeg

    rm -Rf ./build
    mkdir build
    cd build

    # 6) Making VISUAL project with cmake 2.6 - 2.8
    ${CMAKE} -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DUSE_OMP=ON -G "$PIPOL_WIN_COMPILER$SYSTEMOS" ..

    # 7) Setting environment
    if expr $PIPOL_IMAGE : ".*amd64.*"; then
        "$MSVC_PATH/VC/vcvarsall.bat" amd64
    else
        "$MSVC_PATH/Common7/Tools/vsvars32.bat"
    fi

    # 9) Cleaning the project
    "$VCBUILD" OpenMEEG.sln /Clean "Release|$SYSTEMCOMPILE"

    # 10) Project build
    "$VCBUILD" OpenMEEG.sln /Build "Release|$SYSTEMCOMPILE"

    # # 11) Project build
    # "$VCBUILD" INSTALL.vcproj "Release"

    # 12) CDASH submision
    ${CTEST} -D "$1"

    # 13) To make a Windows Package Installer
    ${CPACK} -G "NSIS"
    ${CPACK} -G "ZIP"
    ${CPACK} -G "TGZ"

    # 14) Save the package
    mkdir -p $PIPOL_HOMEDIR/.pipol/packages/openmeeg-${BRANCH}
    cp OpenMEEG-*.* $PIPOL_HOMEDIR/.pipol/packages/openmeeg-${BRANCH}
}
