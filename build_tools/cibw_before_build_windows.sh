#!/bin/bash

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1
echo "Using project root \"${ROOT}\" on RUNNER_OS=\"${RUNNER_OS}\""
cd $ROOT

# Build the Python bindings on Windows
rm -Rf build
cp -a build_nopython build
cmake -B build -DENABLE_PYTHON=ON .
cmake --build build
cp -av build/wrapping/python/openmeeg/* wrapping/python/openmeeg/
