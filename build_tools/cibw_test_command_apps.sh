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
    ldd ./lib/libOpenMEEG.so
    ldd ./lib/libOpenMEEGMaths.so
    # readelf -d binary-or-library | head -20
elif [[ "${RUNNER_OS}" == "macOS" ]]; then
    tar xzfv $ROOT/installers/OpenMEEG-*-*.tar.gz
    cd OpenMEEG-*
    sudo rm -Rf /usr/local/gfortran  # Intel
    sudo rm -Rf /opt/gfortran-darwin-arm64
    otool -L ./lib/libOpenMEEG.dylib
    otool -L ./lib/libOpenMEEGMaths.dylib
elif [[ "${RUNNER_OS}" == "Windows" ]]; then
    ROOT=$(cygpath -u $ROOT)
    tar xzfv $ROOT/installers/OpenMEEG-*.tar.gz
    cd OpenMEEG-*
    $ROOT/Dependencies/Dependencies.exe -modules $(cygpath -w $PWD/bin/om_minverser.exe)
else
    echo "Unknown RUNNER_OS=\"${RUNNER_OS}\""
    exit 1
fi
./bin/om_minverser --help
./bin/om_assemble --help
