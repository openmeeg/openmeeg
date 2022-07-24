#!/bin/bash

if [[ "$1" != "" ]]; then
    export OPENMEEG_DATA_PATH=$1/data
elif [[ "$OPENMEEG_DATA_PATH" == "" ]] && [[ -f "$PWD/data/Head1/Head1.geom" ]]; then
    export OPENMEEG_DATA_PATH=$PWD/data
fi
# TODO: We need to reenable the tests that this comments out. Something about
# how cmake SWIG wraps is better than how setuptools does it
if [[ "$OPENMEEG_BAD_TYPE" == "" ]]; then
    export OPENMEEG_BAD_TYPE=1
fi
if [[ "$OPENMEEG_BAD_MSVC" == "" ]]; then
    export OPENMEEG_BAD_MSVC=$(python -c "import sys; print(int(sys.platform=='win32'))")
fi
TEST_PATH=$(python -c 'from pathlib import Path; import openmeeg; print(Path(openmeeg.__file__).parent)')
echo "OPENMEEG_DATA_PATH=\"$OPENMEEG_DATA_PATH\""
echo "OPENMEEG_BAD_TYPE=\"$OPENMEEG_BAD_TYPE\""
echo "OPENMEEG_BAD_MSVC=\"$OPENMEEG_BAD_MSVC\""
echo "OPENMEEG_BAD_PYPY=\"$OPENMEEG_BAD_PYPY\""
echo
set -xe
pytest -rav "$TEST_PATH"
