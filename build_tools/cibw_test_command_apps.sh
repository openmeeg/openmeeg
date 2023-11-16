#!/bin/bash

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
ROOT=$1

# For now we don't do this because the dependencies are not bundled
# at this stage
# pushd $ROOT/build
# ctest -C release || ctest -C release --rerun-failed --output-on-failure
# popd

ls -al $ROOT
echo
echo "Installers:"
ls -al $ROOT/installers
echo

if [[ "${RUNNER_OS}" == "Linux" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    # readelf -d binary-or-library | head -20
    ./bin/om_minverser --help
    # ldd -v ./lib/libOpenMEEGMaths.so.1
    # ldd -v ./lib/libOpenMEEG.so
elif [[ "${RUNNER_OS}" == "macOS" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    sudo rm /usr/local/gfortran/lib/*.dylib
    otool -L ./lib/libgfortran.*.dylib
    otool -L ./lib/libOpenMEEGMaths.dylib
    otool -L ./lib/libOpenMEEG.dylib
    ./bin/om_minverser --help
elif [[ "${RUNNER_OS}" == "Windows" ]]; then
    ROOT=$(cygpath -u $ROOT)
    tar xzfv $ROOT/installers/OpenMEEG-*.tar.gz
    cd OpenMEEG-*
    $ROOT/Dependencies/Dependencies.exe -modules $(cygpath -w $PWD/bin/om_minverser.exe)
    ./bin/om_minverser --help
else
    echo "Unknown RUNNER_OS=\"${RUNNER_OS}\""
    exit 1
fi
