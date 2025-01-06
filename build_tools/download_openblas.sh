#!/bin/bash -ef

BLAS_VER="v0.3.20-140-gbfd9c1b5"
if [[ "$1" == "" ]]; then
    PLATFORM=$(python -c "import sys; print(sys.platform)")
else
    PLATFORM=$1
fi

# https://anaconda.org/multibuild-wheels-staging/openblas-libs/files
if [[ "${PLATFORM}" == 'win'* ]]; then
    BLAS_FILENAME="openblas-${BLAS_VER}-win_amd64-gcc_8_1_0.zip"
elif [[ "${PLATFORM}" == 'linux'* ]]; then
    BLAS_FILENAME="openblas-${BLAS_VER}-manylinux2014_x86_64.tar.gz"
else
    echo "Unknown/unsupported PLATFORM=\"${PLATFORM}\""
    exit 1
fi

echo "Downloading OpenBLAS and setting cmake flags for PLATFORM=\"${PLATFORM}\""
curl -LO https://anaconda.org/multibuild-wheels-staging/openblas-libs/${BLAS_VER}/download/${BLAS_FILENAME}
if [[ "${PLATFORM}" == 'win'* ]]; then
    unzip ${BLAS_FILENAME} -fod openblas
    export OPENBLAS_LIB=${PWD}/openblas/64/lib
    export OPENBLAS_INCLUDE=${PWD}/openblas/64/include
    BLAS_EXT="${BLAS_VER}-gcc_8_1_0"
    pushd $OPENBLAS_LIB
    # TODO: This conditional should work but it does not...
    # if [[ "${{ matrix.blas_linking }}" != "static" ]]; then
    cp ../bin/libopenblas_${BLAS_EXT}.dll .
    cp libopenblas_${BLAS_EXT}.dll openblas.dll
    mv libopenblas_${BLAS_EXT}.dll.a openblas.dll.a
    # fi
    mv libopenblas_${BLAS_EXT}.a openblas.a
    mv libopenblas_${BLAS_EXT}.def openblas.def
    mv libopenblas_${BLAS_EXT}.exp openblas.exp
    mv libopenblas_${BLAS_EXT}.lib openblas.lib
    popd
    export LIB="$(cygpath -w $OPENBLAS_LIB);${LIB}"
    export CMAKE_CXX_FLAGS="-I$(cygpath -m $OPENBLAS_INCLUDE)"
elif [[ "${PLATFORM}" == 'linux'* ]]; then
    mkdir openblas
    tar xzfv ${BLAS_FILENAME} -C openblas
    BLAS_DIR=${PWD}/openblas/usr/local
    export OPENBLAS_LIB=${BLAS_DIR}/lib
    export OPENBLAS_INCLUDE=${PWD}/openblas/usr/local/include
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE -L$OPENBLAS_LIB"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
else
    echo "Unknown PLATFORM=${PLATFORM}"
    exit 1
fi
