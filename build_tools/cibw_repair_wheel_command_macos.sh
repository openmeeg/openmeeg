#!/bin/bash

set -e
if [[ "$1" == "" ]] || [[ "$2" == "" ]] || [[ "$3" == "" ]]; then
    echo "Usage: $0 <DELOCATE_ARCHS> <DEST_DIR> <WHEEL>"
    exit 1
fi
DELOCATE_ARCHS=$1
DEST_DIR=$2
WHEEL=$3
set -x

export DYLD_LIBRARY_PATH="$GITHUB_WORKSPACE/install/lib"
delocate-listdeps "$WHEEL"
delocate-wheel --require-archs "$DELOCATE_ARCHS" -w "$DEST_DIR" "$WHEEL"
