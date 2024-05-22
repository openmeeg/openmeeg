#!/bin/bash

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1
echo "Using project root \"${ROOT}\" on RUNNER_OS=\"${RUNNER_OS}\""
cd $ROOT
pwd

# Build the Python bindings on Windows
ls -al
rm -Rf build
cp -a build_nopython build
which python
python --version
python -m pip install --upgrade --only-binary="numpy" "numpy>=2.0.0rc1,<3" "setuptools>=68.0.0" "setuptools_scm>=6.2" "wheel>=0.37.0"
cmake -B build -DENABLE_PYTHON=ON -DPython3_EXECUTABLE="$(which python)" .
cmake --build build --config Release
cp -av build/wrapping/python/openmeeg/*.pyd build/wrapping/python/openmeeg/_openmeeg_wrapper.py wrapping/python/openmeeg/
rm -Rf build
git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)"
