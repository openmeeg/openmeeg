#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1
echo "Using project root \"${ROOT}\" on RUNNER_OS=\"${RUNNER_OS}\""
cd $ROOT
pwd

# Let's have NumPy help us out, but we need to tell it to build for the correct
# macOS platform
if [[ "$CIBW_ARCHS_MACOS" == "arm64" ]]; then
    export _PYTHON_HOST_PLATFORM="macosx-11.0-arm64"
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

# PLATFORM can be:
# linux-x86_64
# macosx-x86_64
# macosx-arm64
# win-amd64

if [[ "$PLATFORM" == "linux-x86_64" ]]; then
    dnf -y install epel-release
    dnf -y install hdf5-devel matio-devel
    export OPENBLAS_INCLUDE=/usr/local/include
    export OPENBLAS_LIB=/usr/local/lib
    export CMAKE_CXX_FLAGS="-lgfortran -lpthread -I$OPENBLAS_INCLUDE"
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
elif [[ "$PLATFORM" == 'macosx-'* ]]; then
    brew install boost swig
    BLAS_DIR=/usr/local
    OPENBLAS_INCLUDE=$BLAS_DIR/include
    OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE -L$OPENBLAS_LIB"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
    echo "Building for CIBW_ARCHS_MACOS=\"$CIBW_ARCHS_MACOS\""
    if [[ "$CIBW_ARCHS_MACOS" == "x86_64" ]]; then
        export VCPKG_DEFAULT_TRIPLET="x64-osx-release-10.9"
    elif [[ "$CIBW_ARCHS_MACOS" == "arm64" ]]; then
        export VCPKG_DEFAULT_TRIPLET="arm64-osx-release-10.9"
        CMAKE_OSX_ARCH_OPT="-DCMAKE_OSX_ARCHITECTURES=arm64"
    else
        echo "Unknown CIBW_ARCHS_MACOS=\"$CIBW_ARCHS_MACOS\""
        exit 1
    fi
    if [[ "$CIBW_ARCHS_MACOS" == "arm64" ]]; then
        # The deps were compiled locally on 2022/07/19 on an M1 machine and uploaded
        curl -L https://osf.io/download/x45fz > openmeeg-deps-arm64-osx-release-10.9.tar.gz
        tar xzfv openmeeg-deps-arm64-osx-release-10.9.tar.gz
        CMAKE_PREFIX_PATH_OPT="-DCMAKE_PREFIX_PATH=$PWD/vcpkg_installed/arm64-osx-release-10.9"
    else
        source ./build_tools/setup_vcpkg_compilation.sh
    fi
    CMAKE_OSX_ARCH_OPT="-DCMAKE_OSX_ARCHITECTURES=${CIBW_ARCHS_MACOS}"
    OPENMP_OPT="-DUSE_OPENMP=OFF"
elif [[ "$PLATFORM" == "win-amd64" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows-release-static"
    export CMAKE_GENERATOR="Visual Studio 16 2019"
    source ./build_tools/setup_vcpkg_compilation.sh
    source ./build_tools/download_openblas.sh windows  # NumPy doesn't install the headers for Windows
    pip install delvewheel
    SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
    #export BLA_STATIC_OPT="-DBLA_STATIC=ON"
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export DISABLE_CCACHE=1
pip install cmake
./build_tools/cmake_configure.sh -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${ROOT}/install ${OPENMP_OPT} ${SYSTEM_VERSION_OPT} ${CMAKE_OSX_ARCH_OPT} -DENABLE_APPS=OFF ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT} ${LAPACK_LIBRARIES_OPT}
cmake --build build --target install --config release
# make life easier for auditwheel/delocate/delvewheel
if [[ "$PLATFORM" == 'linux'* ]]; then
    ls -al install/lib64/*.so*
    cp install/lib64/*.so* /usr/local/lib/
elif [[ "$PLATFORM" == 'macosx'* ]]; then
    ls -al install/lib/*.dylib*
    sudo mkdir -p /usr/local/lib
    sudo cp install/lib/*.dylib* /usr/local/lib/
else
    ls -al $PWD/install/bin/*.dll*
    cp $PWD/install/bin/*.dll* .
fi
