#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]]; then
    echo "Usage: $0 <DEST_DIR> <WHEEL>"
    exit 1
fi
DEST_DIR=$1
WHEEL=$2
set -x

export PATH="(cygwin -u $GITHUB_WORKSPACE)/install/bin:$PATH"
delvewheel repair -w "$DEST_DIR" "$WHEEL"
