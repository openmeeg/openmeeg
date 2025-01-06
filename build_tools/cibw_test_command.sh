#!/bin/bash

if [[ "$1" != "" ]]; then
    export OPENMEEG_DATA_PATH=$1/data
elif [[ "$OPENMEEG_DATA_PATH" == "" ]] && [[ -f "$PWD/data/Head1/Head1.geom" ]]; then
    export OPENMEEG_DATA_PATH=$PWD/data
fi
# TODO: We need to reenable the tests that this comments out. Something about
# how cmake SWIG wraps is better than how setuptools does it
echo "OPENMEEG_DATA_PATH=\"$OPENMEEG_DATA_PATH\""
echo
set -xe
python -m threadpoolctl -i numpy openmeeg
echo ""
pytest --fixtures "$TEST_PATH"
echo ""
pwd
echo ""
pytest --tb=short -ra -m "not slow" -vv --pyargs openmeeg
echo ""
# Smoke test for https://github.com/swig/swig/issues/3061
PYTHONFAULTHANDLER=1 PYTHONWARNINGS=error python -uc "import openmeeg"
