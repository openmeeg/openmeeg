if [[ "$USE_PROJECT" == "0" ]]; then
    cd OpenMEEG
fi

mkdir build
cd build

set args = ""

if [[ "$USE_PROJECT" == "1" ]]; then
    # Build OpenMEEG and use installed versions of libraries for lapack and matio.
    set args = "-DUSE_SYSTEM_clapack:BOOL=${USE_SYSTEM} \
                -DUSE_SYSTEM_matio:BOOL=${USE_SYSTEM} \
                -DUSE_SYSTEM_hdf5:BOOL=${USE_SYSTEM} \
                -DUSE_SYSTEM_zlib:BOOL=${USE_SYSTEM}"
fi

cmake \
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
    ${args} \
    -DCMAKE_SKIP_RPATH:BOOL=OFF \
    ..

#cat CMakeCache.txt

make
