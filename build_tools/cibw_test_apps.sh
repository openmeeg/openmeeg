#!/bin/bash

# Check that an app installer tarball is self-contained: run on a machine (or,
# on Linux, in a bare container) that has none of the build dependencies, so
# any library the tarball fails to bundle cannot be picked up from elsewhere.

set -eo pipefail
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <TARBALL>"
    exit 1
fi
TARBALL=$(realpath "$1")
WORK=$(mktemp -d)
tar xzf "$TARBALL" -C "$WORK"
cd "$WORK"/OpenMEEG-*
ls -alR bin lib

BAD=""
case "$(uname -s)" in
    Linux)
        for f in bin/* lib/*.so*; do
            MISSING=$(ldd "$f" | grep "not found" || true)
            if [[ "$MISSING" != "" ]]; then
                BAD="${BAD}${f}:"$'\n'"${MISSING}"$'\n'
            fi
        done
        ;;
    Darwin)
        # Anything but system libraries must be found relative to the bundle
        for f in bin/* lib/*.dylib; do
            OUTSIDE=$(otool -L "$f" | tail -n +2 | awk '{print $1}' | grep -Ev '^(/usr/lib/|/System/|@rpath/|@loader_path/|@executable_path/)' || true)
            if [[ "$OUTSIDE" != "" ]]; then
                BAD="${BAD}${f}:"$'\n'"${OUTSIDE}"$'\n'
            fi
        done
        ;;
esac
if [[ "$BAD" != "" ]]; then
    echo "Unbundled libraries:"
    echo "$BAD"
    exit 1
fi

# Actually loading them is the real test (on Windows, the only one)
./bin/om_minverser --help
./bin/om_assemble --help
