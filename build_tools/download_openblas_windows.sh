#!/bin/bash -ef

curl -LO https://anaconda.org/multibuild-wheels-staging/openblas-libs/v0.3.18/download/openblas-v0.3.18-win_amd64-gcc_8_1_0.zip
BLAS_VER=v0.3.18
BLAS_EXT=${BLAS_VER}-gcc_8_1_0
unzip openblas-${BLAS_VER}-win_amd64-gcc_8_1_0.zip -d openblas
export OPENBLAS_LIB=${PWD}/openblas/64/lib
export OPENBLAS_INCLUDE=${PWD}/openblas/64/include
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
