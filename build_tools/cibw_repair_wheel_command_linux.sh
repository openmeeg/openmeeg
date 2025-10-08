#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]]; then
    echo "Usage: $0 <DEST_DIR> <WHEEL>"
    exit 1
fi
DEST_DIR=$1
WHEEL=$2
set -x

export LD_LIBRARY_PATH="/project/install/lib:$LD_LIBRARY_PATH"
auditwheel show "$WHEEL"
auditwheel repair -w "$DEST_DIR" "$WHEEL"
