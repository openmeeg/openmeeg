#!/bin/bash -ef

which gcc
gcc --version
INSTALL_PREFIX=/c/opt/OpenBLAS
export OPENBLAS_LIB=$INSTALL_PREFIX/lib
export OPENBLAS_INCLUDE=$INSTALL_PREFIX/include
if [ -f "${OPENBLAS_LIB}/openblas.a" ]; then
    echo "Completed OpenBLAS build found, returning..."
    return  # exit as a sourced script without killing the calling shell
fi
if [ -z "$DYNAMIC_ARCH" ]; then
    DYNAMIC_ARCH=1
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
# TODO: This should be cross-checked with
# https://github.com/MacPython/openblas-libs/blob/master/tools/build_openblas.sh#L69-L74
make BINARY=32 DYNAMIC_ARCH=$DYNAMIC_ARCH USE_THREAD=1 USE_OPENMP=0 \
     NUM_THREADS=24 NO_WARMUP=1 NO_AFFINITY=1 CONSISTENT_FPCSR=1 \
     BUILD_LAPACK_DEPRECATED=1 TARGET=PRESCOTT BUFFERSIZE=20 \
     COMMON_OPT="$cflags" \
     FCOMMON_OPT="$fflags" \
     MAX_STACK_ALLOC=2048
make install PREFIX=$INSTALL_PREFIX
test -d $OPENBLAS_LIB
ls -al $OPENBLAS_LIB
test -f $OPENBLAS_LIB/libopenblas.a
cp $OPENBLAS_LIB/libopenblas.a $OPENBLAS_LIB/openblas.a
cp $OPENBLAS_LIB/libopenblas.dll $OPENBLAS_LIB/openblas.dll
test -d $OPENBLAS_INCLUDE
ls -al $OPENBLAS_INCLUDE
test -f $OPENBLAS_INCLUDE/cblas.h
