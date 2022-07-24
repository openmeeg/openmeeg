#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]]; then
    echo "Usage: $0 <DEST_DIR> <WHEEL>"
    exit 1
fi
DEST_DIR=$1
WHEEL=$2

ls -al $GITHUB_WORKSPACE/install/bin
set -x
delvewheel repair --add-path="$(cygpath -w $GITHUB_WORKSPACE/install/bin)" -w "$DEST_DIR" "$WHEEL"
