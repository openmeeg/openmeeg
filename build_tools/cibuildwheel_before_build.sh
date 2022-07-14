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
ls -al .
ls -al build_tools/
if [[ "$PLATFORM" == "linux" ]]; then
    apt-get -y install liblapacke-dev libhdf5-dev libmatio-dev libopenblas-dev
elif [[ "$PLATFORM" == "darwin" ]]; then
    brew install hdf5 libmatio boost swig openblas
    BLAS_DIR=/usr/local/opt/openblas
    OPENBLAS_INCLUDE=$BLAS_DIR/include
    OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE -L$OPENBLAS_LIB"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
elif [[ "$PLATFORM" == "win32" ]]; then
    source ./build_tools/setup_windows_compilation.sh
    source ./build_tools/download_openblas_windows.sh
    pip install delvewheel
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT=-DENABLE_PYTHON=OFF
export STATIC_OPT=-DBLA_STATIC=ON
pip install cmake
./build_tools/cmake_configure.sh -DCMAKE_INSTALL_PREFIX=${ROOT}/install
cmake --build build --target install
