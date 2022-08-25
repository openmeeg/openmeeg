#!/bin/bash

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1

set -e
ls -al $ROOT
echo
echo "Installers:"
ls -alR $ROOT/installers
echo

if [[ "${RUNNER_OS}" == "Linux" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    ldd ./lib/libOpenMEEG.so
    ldd ./lib/libOpenMEEGMaths.so
    # readelf -d binary-or-library | head -20
    ./bin/om_minverser --help
elif [[ "${RUNNER_OS}" == "macOS" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    otool -L ./lib/libOpenMEEG.dylib
    otool -L ./lib/libOpenMEEGMaths.dylib
    ./bin/om_minverser --help
elif [[ "${RUNNER_OS}" == "Windows" ]]; then
    ROOT=$(cygpath -u $ROOT)
    tar xzfv $ROOT/installers/OpenMEEG-*.tar.gz
    cd OpenMEEG-*
    ./bin/om_minverser --help
else
    echo "Unknown RUNNER_OS=\"${RUNNER_OS}\""
    exit 1
fi
