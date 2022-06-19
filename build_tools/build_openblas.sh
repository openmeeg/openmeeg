#!/bin/bash -ef

which gcc
gcc --version
INSTALL_PREFIX=/c/opt/OpenBLAS
export OPENBLAS_LIB=$INSTALL_PREFIX/lib
export OPENBLAS_INCLUDE=$INSTALL_PREFIX/include
if ! -f "${OPENBLAS_LIB}/openblas.a"; then
    echo "Completed OpenBLAS build found, returning..."
    exit 0
fi
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
git submodule update --init --recursive
git fetch --all --tags --prune
git checkout -b v0.3.20 tags/v0.3.20
git clean -fxd
git reset --hard
cflags="-O2 -march=x86-64 -mtune=generic -fno-asynchronous-unwind-tables"
fflags="$cflags -frecursive -ffpe-summary=invalid,zero"
make BINARY=32 DYNAMIC_ARCH=1 USE_THREAD=1 USE_OPENMP=0 \
    NO_WARMUP=1 BUILD_LAPACK_DEPRECATED=1 \
    COMMON_OPT="$cflags" FCOMMON_OPT="$fflags"
make install PREFIX=$INSTALL_PREFIX
test -d $OPENBLAS_LIB
ls -al $OPENBLAS_LIB
test -f $OPENBLAS_LIB/libopenblas.a
cp $OPENBLAS_LIB/libopenblas.a $OPENBLAS_LIB/openblas.a
cp $OPENBLAS_LIB/libopenblas.dll $OPENBLAS_LIB/openblas.dll
test -d $OPENBLAS_INCLUDE
ls -al $OPENBLAS_INCLUDE
test -f $OPENBLAS_INCLUDE/cblas.h
