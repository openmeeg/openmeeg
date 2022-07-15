#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1
PLATFORM=$(python -c "import sys; print(sys.platform)")
echo "Using project root for platform \"${PLATFORM}\": \"${ROOT}\""
cd $ROOT
pwd
# TODO: Use newer OpenBLAS on Linux by downloading it!
if [[ "$PLATFORM" == "linux" ]]; then
    apt-get update -q
    apt-get -yq install libhdf5-dev libmatio-dev libboost-dev
    source ./build_tools/download_openblas.sh linux
    LAPACK_LIBRARIES_OPT="-DLAPACK_LIBRARIES=$OPENBLAS_LIB/libopenblas.a"
elif [[ "$PLATFORM" == "darwin" ]]; then
    brew install hdf5 libmatio boost swig openblas
    BLAS_DIR=/usr/local/opt/openblas
    OPENBLAS_INCLUDE=$BLAS_DIR/include
    OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE -L$OPENBLAS_LIB"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
elif [[ "$PLATFORM" == "win32" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows"
    export CMAKE_GENERATOR="Visual Studio 16 2019"
    source ./build_tools/setup_windows_compilation.sh
    source ./build_tools/download_openblas.sh windows
    pip install delvewheel
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT=-DENABLE_PYTHON=OFF
export BLA_STATIC_OPT=-DBLA_STATIC=ON
export BLA_IMPLEMENTATION=OpenBLAS
export DISABLE_CCACHE=1
pip install cmake
./build_tools/cmake_configure.sh -DCMAKE_INSTALL_PREFIX=${ROOT}/install -DVCPKG_BUILD_TYPE=release -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${LAPACK_LIBRARIES_OPT}
cmake --build build --target install
