#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]]; then
    echo "Usage: $0 <DEST_DIR> <WHEEL>"
    exit 1
fi
DEST_DIR=$1
WHEEL=$2

ADD_PATH=$GITHUB_WORKSPACE/install/bin
ls -al $GITHUB_WORKSPACE/install/bin
ADD_PATH=$(cygpath -w $ADD_PATH)
echo "Adding path \"$ADD_PATH\""
set -x
delvewheel show --add-path="$ADD_PATH" "$WHEEL"
delvewheel repair --add-path="$ADD_PATH" -w "$DEST_DIR" "$WHEEL"
