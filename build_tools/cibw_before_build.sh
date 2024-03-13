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

IS_PYPY=$(python -c "import sys; print(int('pypy' in sys.implementation.name))")
python -m pip install --upgrade --pre --only-binary="numpy" "numpy>=2.0.0b1,<3"
pip install "setuptools>=68.0.0" "wheel>=0.37.0"
