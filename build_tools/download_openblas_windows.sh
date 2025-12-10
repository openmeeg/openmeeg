#!/bin/bash -ef

set -exo pipefail

BLAS_FILENAME="OpenBLAS-0.3.30-x64.zip"
curl -LO https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.30/${BLAS_FILENAME}
unzip ${BLAS_FILENAME} -d openblas
ls -alR openblas
export OPENBLAS_LIB=${PWD}/openblas/lib
export OPENBLAS_INCLUDE=${PWD}/openblas/include
pushd $OPENBLAS_LIB
cp ../bin/libopenblas.dll .
cp libopenblas.dll openblas.dll
mv libopenblas.dll.a openblas.dll.a
mv libopenblas.a openblas.a
mv libopenblas.def openblas.def
mv libopenblas.exp openblas.exp
mv libopenblas.lib openblas.lib
popd
export LIB="$(cygpath -w $OPENBLAS_LIB);${LIB}"
export CMAKE_CXX_FLAGS="-I$(cygpath -m $OPENBLAS_INCLUDE)"
