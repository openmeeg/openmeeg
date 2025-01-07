#!/bin/bash

if [[ "$1" != "" ]]; then
    export OPENMEEG_DATA_PATH=$1/data
elif [[ "$OPENMEEG_DATA_PATH" == "" ]] && [[ -f "$PWD/data/Head1/Head1.geom" ]]; then
    export OPENMEEG_DATA_PATH=$PWD/data
fi
echo "OPENMEEG_DATA_PATH=\"$OPENMEEG_DATA_PATH\""
set -xe
python -m threadpoolctl -i numpy openmeeg
pwd
pytest --tb=short -ra -m "not slow" -vv --pyargs openmeeg
# Smoke test for https://github.com/swig/swig/issues/3061
PYTHONFAULTHANDLER=1 PYTHONWARNINGS=error python -uc "import openmeeg"

# Rerun the "bad" way
TEST_PATH=$(python -c 'from pathlib import Path; import openmeeg; print(Path(openmeeg.__file__).parent)')￼
pytest --tb=short -ra -m "not slow" -vv "$TEST_PATH"
