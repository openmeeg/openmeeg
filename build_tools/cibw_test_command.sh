#!/bin/bash

if [[ "$1" != "" ]]; then
    export OPENMEEG_DATA_PATH=$1/data
elif [[ "$OPENMEEG_DATA_PATH" == "" ]] && [[ -f "$PWD/data/Head1/Head1.geom" ]]; then
    export OPENMEEG_DATA_PATH=$PWD/data
fi
# TODO: We need to reenable the tests that this comments out. Something about
# how cmake SWIG wraps is better than how setuptools does it
TEST_PATH=$(python -c 'from pathlib import Path; import openmeeg; print(Path(openmeeg.__file__).parent)')
echo "OPENMEEG_DATA_PATH=\"$OPENMEEG_DATA_PATH\""
echo
set -xe
python -m threadpoolctl -i numpy openmeeg
pytest --tb=short -ra -m "not slow" -vv "$TEST_PATH"
