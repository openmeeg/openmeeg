#!/bin/bash -ef

# Go to the repo root
DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
cd $DIR/..

if [[ "$VCPKG_DEFAULT_TRIPLET" == "" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows"
fi

if [[ "$VCPKG_DEFAULT_TRIPLET" == "x64-mingw-dynamic" ]]; then
    export CMAKE_GENERATOR="MinGW Makefiles"
    export LINKER_OPT="-s"
elif [[ "$VCPKG_DEFAULT_TRIPLET" == "x64-windows" ]]; then
    if [[ "$CMAKE_GENERATOR" == "" ]]; then  # assume we're using an old version
        CMAKE_GENERATOR="Visual Studio 15 2017"
    fi
    export CMAKE_GENERATOR="$CMAKE_GENERATOR"
    export CMAKE_GENERATOR_PLATFORM="x64"
    export SDK_OPT="-DCMAKE_SYSTEM_VERSION=8.1"
    export TOOLSET_OPT="-DCMAKE_GENERATOR_TOOLSET=v141"
else
    echo "Unknown VCPKG_DEFAULT_TRIPLET: '${VCPKG_DEFAULT_TRIPLET}'"
    exit 1
fi

if [ ! -d vcpkg ]; then
    echo "Getting vcpkg..."
    git clone https://github.com/Microsoft/vcpkg.git --depth=1
    cd vcpkg
    git fetch origin 2022.05.10:use --depth=1
    git checkout use
    ./bootstrap-vcpkg.sh
    cd ..
fi
export VCPKG_INSTALLED_DIR=$(cygpath -m "${PWD}/build/vcpkg_installed")
export VCPKG_INSTALL_OPTIONS="--x-install-root=$VCPKG_INSTALLED_DIR --triplet=$VCPKG_DEFAULT_TRIPLET"
export CMAKE_TOOLCHAIN_FILE=$(cygpath -m "${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake")
if [[ "${PYTHON_OPT}" == "" ]]; then
    export PYTHON_OPT="-DENABLE_PYTHON=ON"
    export PYTHON_INFO_OPT="-DPython3_EXECUTABLE=$(which python)"
fi

if [[ $GITHUB_ENV != "" ]]; then
    echo "VCPKG_INSTALLED_DIR=$VCPKG_INSTALLED_DIR" >> $GITHUB_ENV
    echo "VCPKG_DEFAULT_TRIPLET=$VCPKG_DEFAULT_TRIPLET" >> $GITHUB_ENV
    echo "VCPKG_DEFAULT_HOST_TRIPLET=$VCPKG_DEFAULT_TRIPLET" >> $GITHUB_ENV
    echo "VCPKG_INSTALL_OPTIONS=$VCPKG_INSTALL_OPTIONS" >> $GITHUB_ENV
    echo "CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}" >> $GITHUB_ENV
    echo "SDK_OPT=${SDK_OPT}" >> $GITHUB_ENV
    echo "STRIP_OPT=${STRIP_OPT}" >> $GITHUB_ENV
    echo "TOOLSET_OPT=${TOOLSET_OPT}" >> $GITHUB_ENV
    echo "LINKER_OPT=${LINKER_OPT}" >> $GITHUB_ENV
    echo "CMAKE_GENERATOR=${CMAKE_GENERATOR}" >> $GITHUB_ENV
    echo "CMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}" >> $GITHUB_ENV
fi
test -f $(cygpath -u "$CMAKE_TOOLCHAIN_FILE")
