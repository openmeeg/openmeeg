#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]]; then
    echo "Usage: $0 <DEST_DIR> <WHEEL>"
    exit 1
fi
DEST_DIR=$1
WHEEL=$2

set -x
ADD_PATH=/project/install/lib64
ls -alR "$ADD_PATH"
export LD_LIBRARY_PATH="$ADD_PATH:$LD_LIBRARY_PATH"
auditwheel show "$WHEEL"
auditwheel repair -w "$DEST_DIR" "$WHEEL"
