#!/bin/bash

# 2) Creation of temporary repertory
mkdir /cygdrive/c/Temp

# 3) Configuration of the variables for matlab
export TEMP=/cygdrive/c/Temp
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

# 4) Detectioon of cmake version
if [ -e "/cygdrive/c/CMake 2.8/bin/cmake.exe" ]; then
    VERSION="2.8"
else
    VERSION="2.6"
fi

# 5) Cleaning the projet
cd ${TEMP}

svn checkout svn://scm.gforge.inria.fr/svn/openmeeg/trunk openmeeg-trunk
svn checkout svn://scm.gforge.inria.fr/svn/openmeeg/win32addons win32addons

cd openmeeg-trunk
rm -Rf ./build
mkdir build
cd build

# 6) Making VISUAL project with cmake 2.6 - 2.8
/cygdrive/c/CMake\ $VERSION/bin/cmake.exe -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -G "$PIPOL_WIN_COMPILER$SYSTEMOS" ..

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
/cygdrive/c/CMake\ $VERSION/bin/ctest.exe -D "Nightly"

# 13) To make a Windows Package Installer
/cygdrive/c/CMake\ $VERSION/bin/cpack.exe -G "NSIS"

# 14) Save the package
cp OpenMEEG-*.exe $PIPOL_HOMEDIR/packages/
