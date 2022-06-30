#!/bin/bash -ef

if [[ "${BUILD_TYPE}" == "" ]]; then
    BUILD_TYPE=Release
fi

set -x
cmake -B build \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      $BLA_STATIC_OPT \
      $BLA_IMPL \
      $PYTHON_OPT \
      $DOC_OPT \
      $SDK_OPT \
      $TOOLSET_OPT \
      $PYTHON_INFO_OPT \
      $WERROR_OPT \
      -DCMAKE_EXE_LINKER_FLAGS="$LINKER_OPT" \
      -DCMAKE_SHARED_LINKER_FLAGS="$LINKER_OPT" \
      -DCMAKE_MODULE_LINKER_FLAGS="$LINKER_OPT" \
      -DCMAKE_CXX_COMPILER_LAUNCHER="ccache" \
      -DCMAKE_C_COMPILER_LAUNCHER="ccache" \
      -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
      -DCMAKE_TOOLCHAIN_FILE="$CMAKE_TOOLCHAIN_FILE" \
      -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
      -DTEST_HEAD3=ON \
      "$@"  # any additional args to this script
