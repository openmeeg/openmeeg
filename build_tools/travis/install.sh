if [[ "$USE_PROJECT" == "0" ]]; then
    cd OpenMEEG
fi

mkdir build
cd build

args=""

if [[ "$USE_PROJECT" == "1" ]]; then
    # Build OpenMEEG and use installed versions of libraries for lapack and matio.
    args="-DUSE_SYSTEM_clapack:BOOL=${USE_SYSTEM} \
          -DUSE_SYSTEM_matio:BOOL=${USE_SYSTEM} \
          -DUSE_SYSTEM_hdf5:BOOL=${USE_SYSTEM} \
          -DUSE_SYSTEM_zlib:BOOL=${USE_SYSTEM}"
fi

if [[ "$USE_PYTHON" != "1" && "$USE_PROJECT" == "1" ]]; then  # python requires to use shared libraries
    args="${args} \
          -DBUILD_SHARED_LIBS:BOOL=OFF \
          -DBUILD_SHARED_LIBS_OpenMEEG:BOOL=OFF \
          -DBUILD_SHARED_LIBS_hdf5:BOOL=OFF \
          -DBUILD_SHARED_LIBS_matio:BOOL=OFF \
          -DBUILD_SHARED_LIBS_zlib:BOOL=OFF \
          -DCMAKE_SKIP_RPATH:BOOL=OFF"
fi

echo "cmake args: ${args}"

cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION} \
    -DBUILD_TESTING:BOOL=ON \
    -DENABLE_PYTHON:BOOL=${USE_PYTHON} \
    -DENABLE_COVERAGE:BOOL=${USE_COVERAGE} \
    -DENABLE_PACKAGING:BOOL=${ENABLE_PACKAGING} \
    -DUSE_OMP:BOOL=${USE_OMP} \
    -DUSE_VTK:BOOL=${USE_VTK} \
    -DUSE_GIFTI:BOOL=${USE_GIFTI} \
    -DUSE_CGAL:BOOL=${USE_CGAL} \
    -DAPPLE_STANDALONE:BOOL=${APPLE_STANDALONE} \
    -DBLASLAPACK_IMPLEMENTATION:BOOL=${BLASLAPACK_IMPLEMENTATION} \
    -DCMAKE_SKIP_RPATH:BOOL=OFF \
    ${args} \
    ..

make ${MAKEFLAGS}
