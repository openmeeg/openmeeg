#!/bin/bash

if [[ "$1" != "" ]]; then
    export OPENMEEG_DATA_PATH=$1/data
fi
set -xe
pytest $(python -c 'from pathlib import Path; import openmeeg; print(Path(openmeeg.__file__).parent)')