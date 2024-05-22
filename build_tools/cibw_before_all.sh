#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately

set -eo pipefail
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
cd $1
ROOT=$(pwd)
echo "Using project root \"${ROOT}\" on RUNNER_OS=\"${RUNNER_OS}\""

# Let's have NumPy help us out
if [[ "$RUNNER_OS" == "macOS" ]] && [[ $(uname -m) == 'arm64' ]]; then
    echo "Making /usr/local/lib for macOS arm64"
    sudo mkdir -p /usr/local/lib /usr/local/include
    sudo chown $USER /usr/local/lib /usr/local/include
fi
curl -L https://github.com/numpy/numpy/archive/refs/tags/v1.23.1.tar.gz | tar xz numpy-1.23.1
mv numpy-1.23.1/tools .
mv numpy-1.23.1/numpy .  # on Windows, _distributor_init gets modified
echo "Running NumPy tools/wheels/cibw_before_build.sh $1"
chmod +x ./tools/wheels/cibw_before_build.sh
./tools/wheels/cibw_before_build.sh $1
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")
rm -Rf numpy numpy-1.23.1 tools
echo "Using NumPy PLATFORM=\"${PLATFORM}\""

git checkout LICENSE.txt  # This file is modified by NumPy
git status
git clean

# PLATFORM can be:
# linux-x86_64
# macosx-x86_64
# macosx-arm64
# win-amd64

if [[ "$PLATFORM" == "linux-x86_64" ]]; then
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux
    dnf -y install epel-release
    dnf -y install hdf5-devel matio-devel
    export OPENBLAS_INCLUDE=/usr/local/include
    export OPENBLAS_LIB=/usr/local/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export LINKER_OPT="-lgfortran -lpthread"
    export DISABLE_CCACHE=1
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
elif [[ "$PLATFORM" == 'macosx-'* ]]; then
    BLAS_DIR=/usr/local
    OPENBLAS_INCLUDE=$BLAS_DIR/include
    OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
    export LINKER_OPT="-L$OPENBLAS_LIB"
    if [[ "$PLATFORM" == "macosx-x86_64" ]]; then
        export VCPKG_DEFAULT_TRIPLET="x64-osx-release-10.15"
        export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=10.15"
    elif [[ "$PLATFORM" == "macosx-arm64" ]]; then
        export VCPKG_DEFAULT_TRIPLET="arm64-osx-release-11.0"
        export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0"
    else
        echo "Unknown PLATFORM=\"$PLATFORM\""
        exit 1
    fi
    source ./build_tools/setup_vcpkg_compilation.sh
    # libomp can cause segfaults on macos... maybe from version conflicts with OpenBLAS, or from being too recent?
    export OPENMP_OPT="-DUSE_OPENMP=OFF"
    # need SWIG for Python bindings
    brew install swig
elif [[ "$PLATFORM" == "win-amd64" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows-release-static"
    export CMAKE_GENERATOR="Visual Studio 16 2019"
    source ./build_tools/setup_vcpkg_compilation.sh
    source ./build_tools/download_openblas.sh windows  # NumPy doesn't install the headers for Windows
    pip install delvewheel
    export SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export WERROR_OPT="-DENABLE_WERROR=ON"
pip install cmake
./build_tools/cmake_configure.sh -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${ROOT}/install ${CMAKE_PREFIX_PATH_OPT} -DENABLE_APPS=OFF ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT} ${LAPACK_LIBRARIES_OPT}
cmake --build build --target install --target package --config release

# Put DLLs where they can be found
if [[ "$PLATFORM" == 'linux'* ]]; then
    ls -al install/lib64/*.so*
    cp -av install/lib64/*.so* /usr/local/lib/
elif [[ "$PLATFORM" == 'macosx-arm64' ]]; then
    # https://matthew-brett.github.io/docosx/mac_runtime_link.html
    #cp -av $ROOT/vcpkg_installed/arm64-osx-release-11.0/lib/libomp* $ROOT/install/lib/
    otool -L $ROOT/install/lib/libOpenMEEG.1.1.0.dylib
    # install_name_tool -change "@@HOMEBREW_PREFIX@@/opt/libomp/lib/libomp.dylib" "@loader_path/libomp.dylib" $ROOT/install/lib/libOpenMEEG.1.1.0.dylib
    # otool -L $ROOT/install/lib/libOpenMEEG.1.1.0.dylib
elif [[ "$PLATFORM" == 'win'* ]]; then
    cp -av $OPENBLAS_LIB/libopenblas_v0.3.20-140-gbfd9c1b5-gcc_8_1_0.dll install/bin/
fi

# TODO: This is only necessary because SWIG does not work outside cmake yet,
# and we want this on windows
if [[ "$PLATFORM" == 'win'* ]]; then
    mv build build_nopython
fi
echo "ls -al $PWD:"
ls -al

echo "git status:"
git status

echo "git clean"
git clean
