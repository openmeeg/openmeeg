#!/bin/bash -ef

if [[ "${BUILD_TYPE}" == "" ]]; then
    BUILD_TYPE=Release
fi

if [[ "${PYTHON_OPT}" == "" ]]; then
    PYTHON_OPT="-DENABLE_PYTHON=ON"
    PYTHON_EXECUTABLE_OPT="-DPython3_EXECUTABLE=$(which python)"
fi

# Most of the time we want to use ccache, but lets allow for disabling it
# (e.g., on cibuildwheel)
if [[ "${DISABLE_CCACHE}" != "1" ]]; then
    CXX_COMPILER_LAUNCHER_OPT="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"
    C_COMPILER_LAUNCHER_OPT="-DCMAKE_C_COMPILER_LAUNCHER=ccache"
fi

echo "::group:cmake configure"
set -x
cmake -B build \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      $BLA_STATIC_OPT \
      $BLA_IMPL \
      $PYTHON_OPT \
      $PYTHON_EXECUTABLE_OPT \
      $PYTHON_COPY_RUNTIME_DLLS_OPT \
      $DOC_OPT \
      $VTK_OPT \
      $TOOLSET_OPT \
      $WERROR_OPT \
      $CXX_COMPILER_LAUNCHER_OPT \
      $C_COMPILER_LAUNCHER_OPT \
      $VCPKG_TRIPLET_OPT \
      $SYSTEM_VERSION_OPT \
      $SLOW_OPT \
      $OPENMP_OPT \
      -DCMAKE_EXE_LINKER_FLAGS="$LINKER_OPT" \
      -DCMAKE_SHARED_LINKER_FLAGS="$LINKER_OPT" \
      -DCMAKE_MODULE_LINKER_FLAGS="$LINKER_OPT" \
      -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
      -DCMAKE_TOOLCHAIN_FILE="$CMAKE_TOOLCHAIN_FILE" \
      -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
      -DTEST_HEAD3=ON \
      "$@"  # any additional args to this script
set +x
echo "::endgroup"
