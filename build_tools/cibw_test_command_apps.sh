#!/bin/bash

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1

set -e
ls -al $ROOT
ls -alR $ROOT/installers
ls -al $ROOT/installers/OpenMEEG-*-*.*

if [[ "${RUNNER_OS}" == "Linux" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    ./bin/om_minverser --help
elif [[ "${RUNNER_OS}" == "macOS" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    ./bin/om_minverser --help
elif [[ "${RUNNER_OS}" == "Windows" ]]; then
    $ROOT/installers/OpenMEEG-*-*.exe
    om_minverser --help
else
    echo "Unknown RUNNER_OS=\"${RUNNER_OS}\""
    exit 1
fi
