#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately

set -e
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
cd $1
ROOT=$(pwd)
echo "Using project root \"${ROOT}\" on RUNNER_OS=\"${RUNNER_OS}\""

# Let's have NumPy help us out
if [[ "$RUNNER_OS" == "macOS" ]] && [[ $(uname -m) == 'arm64' ]]; then
    echo "Making /usr/local/lib for macOS arm64"
    sudo mkdir -p /usr/local/lib /usr/local/include
    sudo chown $USER /usr/local/lib /usr/local/include
fi
curl -L https://github.com/numpy/numpy/archive/refs/tags/v1.23.1.tar.gz | tar xz numpy-1.23.1
mv numpy-1.23.1/tools .
mv numpy-1.23.1/numpy .  # on Windows, _distributor_init gets modified
echo "Running NumPy tools/wheels/cibw_before_build.sh $1"
chmod +x ./tools/wheels/cibw_before_build.sh
./tools/wheels/cibw_before_build.sh $1
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")
rm -Rf numpy numpy-1.23.1 tools
echo "Using NumPy PLATFORM=\"${PLATFORM}\""

# PLATFORM can be:
# linux-x86_64
# macosx-x86_64
# macosx-arm64
# win-amd64

if [[ "$PLATFORM" == "linux-x86_64" ]]; then
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux
    dnf -y install curl zip unzip tar
    export OPENBLAS_INCLUDE=/usr/local/include
    export OPENBLAS_LIB=/usr/local/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export LINKER_OPT="-lgfortran -lpthread"
    export DISABLE_CCACHE=1
    export VCPKG_DEFAULT_TRIPLET="x64-linux"
    source ./build_tools/setup_vcpkg_compilation.sh
    LAPACK_LIBRARIES_OPT="-DLAPACK_LIBRARIES=/usr/local/lib/libopenblas.a"
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
    LIBDIR_OPT="-DCMAKE_INSTALL_LIBDIR=lib"
elif [[ "$PLATFORM" == 'macosx-'* ]]; then
    BLAS_DIR=/usr/local
    OPENBLAS_INCLUDE=$BLAS_DIR/include
    OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
    export LINKER_OPT="-L$OPENBLAS_LIB"
    export LINKER_OPT="$LINKER_OPT -lgfortran"
    if [[ "$PLATFORM" == "macosx-x86_64" ]]; then
        export VCPKG_DEFAULT_TRIPLET="x64-osx-release-10.15"
        export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=10.15"
        PACKAGE_ARCH_SUFFIX="_Intel"
        # gfortran
        LINKER_OPT="$LINKER_OPT -L/usr/local/gfortran/lib"
        sudo chmod -R a+w /usr/local/gfortran/lib
        name=/usr/local/gfortran/lib
        install_name_tool -change "${name}/libquadmath.0.dylib" "@rpath/libquadmath.0.dylib" ${name}/libgfortran.3.dylib
        install_name_tool -change "${name}/libgcc_s.1.dylib" "@rpath/libgcc_s.1.dylib" ${name}/libgfortran.3.dylib
        install_name_tool -id "@rpath/libgfortran.3.dylib" ${name}/libgfortran.3.dylib
        otool -L ${name}/libgfortran.3.dylib
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=/usr/local/gfortran/lib/libgfortran.3.dylib;/usr/local/gfortran/lib/libquadmath.0.dylib;/usr/local/gfortran/lib/libgcc_s.1.dylib"
    elif [[ "$PLATFORM" == "macosx-arm64" ]]; then
        export VCPKG_DEFAULT_TRIPLET="arm64-osx-release-11.0"
        export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0"
        PACKAGE_ARCH_SUFFIX="_M1"
    else
        echo "Unknown PLATFORM=\"$PLATFORM\""
        exit 1
    fi
    source ./build_tools/setup_vcpkg_compilation.sh
    # libomp can cause segfaults on macos... maybe from version conflicts with OpenBLAS, or from being too recent?
    export OPENMP_OPT="-DUSE_OPENMP=OFF"
    PACKAGE_ARCH_OPT="-DPACKAGE_ARCH_SUFFIX=$PACKAGE_ARCH_SUFFIX"
elif [[ "$PLATFORM" == "win-amd64" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows-release-static"
    export CMAKE_GENERATOR="Visual Studio 16 2019"
    source ./build_tools/setup_vcpkg_compilation.sh
    source ./build_tools/download_openblas.sh windows  # NumPy doesn't install the headers for Windows
    pip install delvewheel
    export SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
    OPENBLAS_DLL=$(ls $OPENBLAS_LIB/libopenblas*.dll)
    echo "OPENBLAS_DLL=\"${OPENBLAS_DLL}\""
    test -f $OPENBLAS_DLL
    LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=$(cygpath -m ${OPENBLAS_DLL})"
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export WERROR_OPT="-DENABLE_WERROR=ON"
pip install cmake
export BLA_STATIC_OPT="-DBLA_STATIC=ON"
./build_tools/cmake_configure.sh -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${ROOT}/install ${LIBDIR_OPT} ${LIBRARIES_INSTALL_OPT} ${PACKAGE_ARCH_OPT} ${CMAKE_PREFIX_PATH_OPT} -DENABLE_APPS=ON ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT} ${LAPACK_LIBRARIES_OPT}
cmake --build build --config release
if [[ "${PLATFORM}" == 'macosx-x86_64'* ]]; then
    for name in OpenMEEG OpenMEEGMaths; do
        install_name_tool -change "/usr/local/gfortran/lib/libgfortran.3.dylib" "@rpath/libgfortran.3.dylib" ./build/${name}/lib${name}.1.1.0.dylib
    done
fi
cmake --build build --target package --target install --config release
mkdir -p installers
cp -av build/OpenMEEG-*-*.* installers/

# Put DLLs where they can be found
if [[ "$PLATFORM" == 'linux'* ]]; then
    ls /
    mkdir -p /output
    cp -av build/OpenMEEG-*-*.* /output/
elif [[ "$PLATFORM" == 'macosx-'* ]]; then
    otool -L $ROOT/build/OpenMEEG/libOpenMEEG.1.1.0.dylib
    otool -L $ROOT/install/lib/libgfortran.*.dylib
else
    ./build_tools/install_dependency_walker.sh
fi
