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
if [[ "$PLATFORM" == "linux" ]]; then
    apt-get update -q
    apt-get -yq install libhdf5-dev libmatio-dev libboost-dev
    source ./build_tools/download_openblas.sh linux
    BLAS_LIBRARIES_OPT="-DBLAS_LIBRARIES=$OPENBLAS_LIB/libopenblas.a"
    LAPACK_LIBRARIES_OPT="-DLAPACK_LIBRARIES=$OPENBLAS_LIB/libopenblas.a"
    HDF5_LIBRARIES_OPT="-DHDF5_LIBRARIES=/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so"
    export CMAKE_CXX_FLAGS="-lgfortran -I$OPENBLAS_INCLUDE"
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
    VCPKG_BUILD_TYPE_OPT="-DVCPKG_BUILD_TYPE=release"
    SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export DISABLE_CCACHE=1
pip install cmake
./build_tools/cmake_configure.sh -DCMAKE_INSTALL_PREFIX=${ROOT}/install ${VCPKG_BUILD_TYPE_OPT} ${SYSTEM_VERSION_OPT} -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT} ${LAPACK_LIBRARIES_OPT} ${HDF5_LIBRARIES_OPT}
cmake --build build --target install
# make life easier for auditwheel/delocate/delvewheel
if [[ "$PLATFORM" == "linux" ]]; then
    ls -al install/lib/*.so*
    cp install/lib/*.so* /usr/local/lib
elif [[ "$PLATFORM" == "darwin" ]]; then
    ls -al install/lib/*.dylib*
    cp install/lib/*.so* /usr/lib
else
    ls -al $PWD/install/bin/*.dll*
    cp $PWD/install/bin/*.dll* .
fi
