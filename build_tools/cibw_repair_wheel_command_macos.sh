#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]] || [[ "$3" == "" ]]; then
    echo "Usage: $0 <DELOCATE_ARCHS> <DEST_DIR> <WHEEL>"
    exit 1
fi
DELOCATE_ARCHS=$1
DEST_DIR=$2
WHEEL=$3

# If running locally, set GITHUB_WORKSPACE to current dir
if [[ "$GITHUB_WORKSPACE" == "" ]]; then
    GITHUB_WORKSPACE=$(pwd)
fi

set -x
ADD_PATH="$GITHUB_WORKSPACE/install/lib"
ls -alR "$ADD_PATH"
export DYLD_LIBRARY_PATH="$ADD_PATH:$DYLD_LIBRARY_PATH"
delocate-listdeps "$WHEEL"
delocate-wheel --require-archs "$DELOCATE_ARCHS" -w "$DEST_DIR" "$WHEEL"
