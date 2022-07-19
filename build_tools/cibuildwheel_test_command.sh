#!/bin/bash

if [[ "$1" != "" ]]; then
    export OPENMEEG_DATA_PATH=$1/data
fi
# TODO: We need to reenable the tests that this comments out. Something about
# how cmake SWIG wraps is better than how setuptools does it
export OPENMEEG_BAD_TYPE=1
export OPENMEEG_BAD_MSVC=$(python -c "import sys; print(int(sys.platform=='win32'))")
set -xe
pytest $(python -c 'from pathlib import Path; import openmeeg; print(Path(openmeeg.__file__).parent)')