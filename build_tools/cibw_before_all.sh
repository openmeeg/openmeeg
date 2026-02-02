#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately

set -eo pipefail
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
KIND=$2
if [[ "$KIND" == "" ]]; then
    KIND=wheel
fi
if [[ "$KIND" != "wheel" ]] && [[ "$KIND" != "app" ]]; then
    echo "Usage $0 <PROJECT_PATH> <KIND>"
    echo "Got KIND=\"$KIND\", should be \"wheel\" or \"app\""
    exit 1
fi
cd $1
ROOT=$(pwd)
if [ -z "$RUNNER_OS" ]; then
    echo "RUNNER_OS is empty, it must be present in the environment"
    exit 1
fi
echo "Using project root \"$ROOT\" on RUNNER_OS=\"${RUNNER_OS}\" to set up KIND=\"$KIND\""

echo "::group::scipy-openblas32"
PLATFORM=$(python -c "import platform; print(f'{platform.system()}-{platform.machine()}')")
echo "PLATFORM=$PLATFORM"
# PLATFORM can be:
# Linux-x86_64
# Linux-aarch64
# Darwin-x86_64
# Darwin-arm64
# Windows-AMD64

python -m pip install "scipy-openblas32!=0.3.30.0.4,!=0.3.30.0.3"
# fix a bug in the headers! https://github.com/OpenMathLib/OpenBLAS/issues/5493
if [[ "$PLATFORM" == 'Darwin-'* ]]; then
    SED_OPT="-i ''"
else
    SED_OPT="-i"
fi
sed $SED_OPT "s/ LAPACKE_/ scipy_LAPACKE_/g" "$(python -c 'import scipy_openblas32; print(scipy_openblas32.get_include_dir())')/lapacke.h"
OPENBLAS_INCLUDE=$(python -c "import scipy_openblas32; print(scipy_openblas32.get_include_dir())")
echo "OPENBLAS_INCLUDE=\"$OPENBLAS_INCLUDE\""
ls -alR $OPENBLAS_INCLUDE
OPENBLAS_LIB_DIR=$(python -c "import scipy_openblas32; print(scipy_openblas32.get_lib_dir())")  # somewhere like "/absolute/path/to/site-packages/scipy_openblas32/lib"
OPENBLAS_LIB_NAME=$(python -c "import scipy_openblas32; print(scipy_openblas32.get_library())")  # typically, "libscipy_openblas32"
# mkdir -p ./.openblas
# echo "./.openblas/scipy_openblas.pc:"
# echo $(python -c "import pathlib, scipy_openblas32; pathlib.Path('./.openblas/scipy_openblas.pc').write_text(scipy_openblas32.get_pkg_config())")
# export PKG_CONFIG_PATH="$PWD/.openblas"
# echo "PKG_CONFIG_PATH=\"$PKG_CONFIG_PATH\""
# cat $PKG_CONFIG_PATH/scipy_openblas.pc
CMAKE_PREFIX_PATH_OPT="-DCMAKE_PREFIX_PATH=$OPENBLAS_LIB_DIR/cmake/openblas"
# export OpenBLAS_DIR="$OPENBLAS_LIB_DIR/cmake/openblas"
echo "CMAKE_MODULE_PATH=\"$CMAKE_MODULE_PATH\""
echo "OPENBLAS_LIB_DIR=\"$OPENBLAS_LIB_DIR\""
ls -alR $OPENBLAS_LIB_DIR
echo "::endgroup::"

git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"

if [[ "$PLATFORM" == 'Linux-'* ]]; then
    echo "::group::yum"
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux
    yum -y install epel-release
    yum -y install curl zip unzip tar ninja-build wget
    echo "::endgroup::"
    echo "::group::matio"
    set -x
    wget https://github.com/tbeu/matio/releases/download/v1.5.30/matio-1.5.30.tar.gz
    tar xvf matio-1.5.30.tar.gz
    pushd matio-1.5.30
    cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
    cmake --build build --config Release --target install
    popd
    set +x
    echo "::endgroup::"
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export DISABLE_CCACHE=1
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
    if [[ "$KIND" == "app" ]]; then
        if [[ "$PLATFORM" == "Linux-x86_64" ]]; then
            export VCPKG_DEFAULT_TRIPLET="x64-linux"
        elif [[ "$PLATFORM" == "Linux-aarch64" ]]; then
            export VCPKG_DEFAULT_TRIPLET="arm64-linux"
        else
            echo "Unknown PLATFORM=\"$PLATFORM\""
            exit 1
        fi
        source ./build_tools/setup_vcpkg_compilation.sh
        LIBDIR_OPT="-DCMAKE_INSTALL_LIBDIR=lib"
        if [[ "$KIND" == "app" ]]; then
            SOS=($OPENBLAS_LIB_DIR/*.so*)
            OLD_IFS=$IFS
            IFS=';'
            SOS=("${SOS[*]}")
            IFS=$OLD_IFS
            LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=${SOS}"
            echo "LIBRARIES_INSTALL_OPT=\"$LIBRARIES_INSTALL_OPT\""
        fi
    fi
    LIB_OUTPUT_DIR="$ROOT/install/lib64"
elif [[ "$PLATFORM" == 'Darwin-'* ]]; then
    if [[ "$PLATFORM" == "Darwin-x86_64" ]]; then
        VC_NAME="x64"
        MIN_VER="10.15"
    elif [[ "$PLATFORM" == "Darwin-arm64" ]]; then
        VC_NAME="arm64"
        MIN_VER="11.0"
    else
        echo "Unknown PLATFORM=\"$PLATFORM\""
        exit 1
    fi
    export VCPKG_DEFAULT_TRIPLET="${VC_NAME}-osx-release-${MIN_VER//./}"  # remove dots
    export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=${MIN_VER}"
    source ./build_tools/setup_vcpkg_compilation.sh
    # libomp can cause segfaults on macos... maybe from version conflicts with OpenBLAS, or from being too recent?
    export OPENMP_OPT="-DUSE_OPENMP=OFF"
    if [[ "$KIND" == "app" ]]; then
        if [[ "$PLATFORM" == "Darwin-x86_64" ]]; then
            PACKAGE_ARCH_SUFFIX="_Intel"
        elif [[ "$PLATFORM" == "Darwin-arm64" ]]; then
            PACKAGE_ARCH_SUFFIX="_M1"
        else
            echo "Unknown PLATFORM=\"$PLATFORM\""
            exit 1
        fi
        PACKAGE_ARCH_OPT="-DPACKAGE_ARCH_SUFFIX=$PACKAGE_ARCH_SUFFIX"
    fi
    LIB_OUTPUT_DIR="$ROOT/install/lib"
elif [[ "$PLATFORM" == "Windows-AMD64" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows-release-static"
    export CMAKE_GENERATOR="Visual Studio 17 2022"
    source ./build_tools/setup_vcpkg_compilation.sh
    pip install delvewheel "pefile!=2024.8.26"
    export SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
    OPENBLAS_DLL="$OPENBLAS_LIB_DIR\\$OPENBLAS_LIB_NAME.dll"
    OPENBLAS_DLL="$(cygpath -u $OPENBLAS_DLL)"
    echo "OPENBLAS_DLL=${OPENBLAS_DLL}"
    test -f $OPENBLAS_DLL
    if [[ "$KIND" == "app" ]]; then
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=$(cygpath -m ${OPENBLAS_DLL})"
    fi
    LIB_OUTPUT_DIR="$ROOT/install/bin"
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export BLAS_LIBRARIES_OPT="-DUSE_SCIPY_OPENBLAS=ON"
export WERROR_OPT="-DENABLE_WERROR=ON"
echo "::group::pip"
pip install --upgrade cmake "swig>=4.2"
echo "::endgroup::"
if [[ "${KIND}" == "wheel" ]]; then
    ./build_tools/cmake_configure.sh -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=$ROOT/install ${CMAKE_PREFIX_PATH_OPT} -DENABLE_APPS=OFF ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT}
    echo "::group::cmake --build"
    cmake --build build --target install --target package --config release
    echo "::endgroup::"
else
    export BLA_STATIC_OPT="-DBLA_STATIC=ON"
    ./build_tools/cmake_configure.sh -DCMAKE_WARN_DEPRECATED=FALSE -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=$ROOT/install ${LIBDIR_OPT} ${LIBRARIES_INSTALL_OPT} ${PACKAGE_ARCH_OPT} ${CMAKE_PREFIX_PATH_OPT} -DENABLE_APPS=ON ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT}
    echo "::group::cmake --build"
    cmake --build build --config release
    echo "::endgroup::"
    echo "::group::cmake --target package"
    cmake --build build --target package --target install --config release
    mkdir -p installers
    cp -av build/OpenMEEG-*-*.* installers/
    echo "::endgroup::"
fi
cp -av "$OPENBLAS_INCLUDE" "$ROOT/install/include/OpenBLAS"  # for build step, need it somewhere we know where it is!
test -f "$ROOT/install/include/OpenBLAS/cblas.h"

# Put DLLs where they can be found
if [[ "$PLATFORM" == 'Linux-'* ]]; then
    mkdir -p /output
    if [[ "$KIND" == "wheel" ]]; then
        ls -al $LIB_OUTPUT_DIR/*.so*
        # For some reason, copying to LIB_OUTPUT_DIR doesn't work for vendoring
        # into the wheel, so copy it somewhere else it will be found.
        # TODO: Maybe someday we can/should just use the EXTRA_INSTALL_LIBRARIES
        # mechanic for this, too
        cp -av "$ROOT/install/lib64"/*.so* /usr/local/lib/
        cp -av "$OPENBLAS_LIB_DIR"/*.so* /usr/local/lib/
    else
        cp -av build/OpenMEEG-*-*.* /output/
    fi
elif [[ "$PLATFORM" == 'Darwin-'* ]]; then
    if [[ "$KIND" == "wheel" ]]; then
        otool -L $LIB_OUTPUT_DIR/libOpenMEEG.1.1.0.dylib
        cp -av "$OPENBLAS_LIB_DIR"/*.dylib* $LIB_OUTPUT_DIR/
    else
        otool -L $ROOT/build/OpenMEEG/libOpenMEEG.1.1.0.dylib
    fi
elif [[ "$PLATFORM" == 'Windows-'* ]]; then
    if [[ "$KIND" == "wheel" ]]; then
        cp -av "$OPENBLAS_DLL" $LIB_OUTPUT_DIR/
    else
        ./build_tools/install_dependency_walker.sh
    fi
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi

if [[ "$KIND" == "wheel" ]]; then
    echo "ls -al $ROOT:"
    ls -al "$ROOT"
fi

echo "git status:"
git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"
