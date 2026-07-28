#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]]; then
    echo "Usage: $0 <DEST_DIR> <WHEEL>"
    exit 1
fi
DEST_DIR=$1
WHEEL=$2

# If running locally, set GITHUB_WORKSPACE to current dir
if [[ "$GITHUB_WORKSPACE" == "" ]]; then
    GITHUB_WORKSPACE=$(pwd)
fi
set -x
ADD_PATH=$GITHUB_WORKSPACE/install/bin
ls -alR "$ADD_PATH"
ADD_PATH=$(cygpath -w $ADD_PATH)
delvewheel show --add-path="$ADD_PATH" "$WHEEL"
delvewheel repair --add-path="$ADD_PATH" --no-mangle="libquadmath-0.dll;libgcc_s_seh-1.dll" -w "$DEST_DIR" "$WHEEL"
