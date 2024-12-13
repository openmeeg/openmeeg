#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately

set -eo pipefail
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi
KIND=$2
if [[ "$KIND" == "" ]]; then
    KIND=wheel
fi
if [[ "$KIND" != "wheel" ]] && [[ "$KIND" != "app" ]]; then
    echo "Usage $0 <PROJECT_PATH> <KIND>"
    echo "Got KIND=\"$KIND\", should be \"wheel\" or \"app\""
    exit 1
fi
cd $1
ROOT=$(pwd)
if [ -z "$RUNNER_OS" ]; then
    echo "RUNNER_OS is empty, it must be present in the environment"
    exit 1
fi
echo "Using project root \"${ROOT}\" on RUNNER_OS=\"${RUNNER_OS}\" to set up KIND=\"$KIND\""

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
git config --global --add safe.directory "$ROOT"
git checkout LICENSE.txt  # This file is modified by NumPy
if [[ "$KIND" == "app" ]]; then
    rm -Rf .ccache gfortran.dmg gfortran-darwin-arm64.tar.gz openblas-v*.zip
fi
git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"

# PLATFORM can be:
# linux-x86_64
# linux-aarch64
# macosx-x86_64
# macosx-arm64
# win-amd64

if [[ "$PLATFORM" == 'linux-'* ]]; then
    echo "::group::yum"
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux
    yum -y install epel-release
    yum -y install hdf5-devel matio-devel curl zip unzip tar ninja-build
    echo "::endgroup::"
    BLAS_DIR=/usr/local
    export OPENBLAS_INCLUDE=$BLAS_DIR/include
    export OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export LINKER_OPT="-lgfortran -lpthread"
    export DISABLE_CCACHE=1
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
    if [[ "$KIND" == "app" ]]; then
        if [[ "$PLATFORM" == "linux-x86_64" ]]; then
            export VCPKG_DEFAULT_TRIPLET="x64-linux"
        elif [[ "$PLATFORM" == "linux-aarch64" ]]; then
            export VCPKG_DEFAULT_TRIPLET="arm64-linux"
        else
            echo "Unknown PLATFORM=\"$PLATFORM\""
            exit 1
        fi
        source ./build_tools/setup_vcpkg_compilation.sh
        LAPACK_LIBRARIES_OPT="-DLAPACK_LIBRARIES=/usr/local/lib/libopenblas.a"
        LIBDIR_OPT="-DCMAKE_INSTALL_LIBDIR=lib"
    fi
elif [[ "$PLATFORM" == 'macosx-'* ]]; then
    BLAS_DIR=/usr/local
    OPENBLAS_INCLUDE=$BLAS_DIR/include
    OPENBLAS_LIB=$BLAS_DIR/lib
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export CMAKE_PREFIX_PATH="$BLAS_DIR"
    export LINKER_OPT="-L$OPENBLAS_LIB"
    if [[ "$PLATFORM" == "macosx-x86_64" ]]; then
        VC_NAME="x64"
        MIN_VER="10.15"
        LIBGFORTRAN="/usr/local/gfortran/lib/libgfortran.3.dylib"
    elif [[ "$PLATFORM" == "macosx-arm64" ]]; then
        VC_NAME="arm64"
        MIN_VER="11.0"
        LIBGFORTRAN="$(find /opt/gfortran-darwin-arm64/lib -name libgfortran.dylib)"
    else
        echo "Unknown PLATFORM=\"$PLATFORM\""
        exit 1
    fi
    export VCPKG_DEFAULT_TRIPLET="${VC_NAME}-osx-release-${MIN_VER}"
    export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=${MIN_VER}"
    source ./build_tools/setup_vcpkg_compilation.sh
    # libomp can cause segfaults on macos... maybe from version conflicts with OpenBLAS, or from being too recent?
    export OPENMP_OPT="-DUSE_OPENMP=OFF"
    if [[ "$KIND" == "app" ]]; then
        GFORTRAN_LIB=$(dirname $LIBGFORTRAN)
        GFORTRAN_NAME=$(basename $LIBGFORTRAN)
        sudo chmod -R a+w $GFORTRAN_LIB
        otool -L $LIBGFORTRAN
        install_name_tool -id "@rpath/${GFORTRAN_NAME}" $LIBGFORTRAN
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=$LIBGFORTRAN"
        if [[ "$PLATFORM" == "macosx-x86_64" ]]; then
            PACKAGE_ARCH_SUFFIX="_Intel"
            install_name_tool -change "${GFORTRAN_LIB}/libgcc_s.1.dylib" "@rpath/libgcc_s.1.dylib" ${LIBGFORTRAN}
            install_name_tool -change "${GFORTRAN_LIB}/libquadmath.0.dylib" "@rpath/libquadmath.0.dylib" ${LIBGFORTRAN}
            LIBRARIES_INSTALL_OPT="$LIBRARIES_INSTALL_OPT;$GFORTRAN_LIB/libgcc_s.1.dylib;$GFORTRAN_LIB/libquadmath.0.dylib"
        else
            PACKAGE_ARCH_SUFFIX="_M1"
            install_name_tool -change "${GFORTRAN_LIB}/libgcc_s.2.dylib" "@rpath/libgcc_s.2.dylib" ${LIBGFORTRAN}
            LIBRARIES_INSTALL_OPT="$LIBRARIES_INSTALL_OPT;$GFORTRAN_LIB/libgcc_s.2.dylib"
            # Doesn't seem like this should be necessary but it is for the arm64 build
            codesign --force -s - $GFORTRAN_LIB/libgcc_s.2.dylib
        fi
        # Need to fix the now-broken signature via ad-hoc signing (at least on arm)
        # https://github.com/matthew-brett/delocate/blob/de38e09acd86b27c795c3d342d132031c45b1aff/delocate/tools.py#L660
        codesign --force -s - $LIBGFORTRAN
        # Set LINKER_OPT after vckpg_compilation.sh because it also sets LINKER_OPT
        export LINKER_OPT="$LINKER_OPT -lgfortran -L$GFORTRAN_LIB"
        PACKAGE_ARCH_OPT="-DPACKAGE_ARCH_SUFFIX=$PACKAGE_ARCH_SUFFIX"
    fi
elif [[ "$PLATFORM" == "win-amd64" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows-release-static"
    export CMAKE_GENERATOR="Visual Studio 17 2022"
    source ./build_tools/setup_vcpkg_compilation.sh
    source ./build_tools/download_openblas.sh windows  # NumPy doesn't install the headers for Windows
    pip install delvewheel "pefile!=2024.8.26"
    export SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
    if [[ "$KIND" == "app" ]]; then
        OPENBLAS_DLL=$(ls $OPENBLAS_LIB/libopenblas*.dll)
        echo "OPENBLAS_DLL=\"${OPENBLAS_DLL}\""
        test -f $OPENBLAS_DLL
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=$(cygpath -m ${OPENBLAS_DLL})"
    fi
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export WERROR_OPT="-DENABLE_WERROR=ON"
echo "::group::pip"
pip install --upgrade cmake "swig>=4.2"
echo "::endgroup::"
if [[ "${KIND}" == "wheel" ]]; then
    ./build_tools/cmake_configure.sh -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${ROOT}/install ${CMAKE_PREFIX_PATH_OPT} -DENABLE_APPS=OFF ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT} ${LAPACK_LIBRARIES_OPT}
    echo "::group::cmake --build"
    cmake --build build --target install --target package --config release
    echo "::endgroup::"
else
    export BLA_STATIC_OPT="-DBLA_STATIC=ON"
    ./build_tools/cmake_configure.sh -DCMAKE_WARN_DEPRECATED=FALSE -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=${ROOT}/install ${LIBDIR_OPT} ${LIBRARIES_INSTALL_OPT} ${PACKAGE_ARCH_OPT} ${CMAKE_PREFIX_PATH_OPT} -DENABLE_APPS=ON ${SHARED_OPT} -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${BLAS_LIBRARIES_OPT} ${LAPACK_LIBRARIES_OPT}
    echo "::group::cmake --build"
    cmake --build build --config release
    if [[ "${PLATFORM}" == 'macosx-'* ]]; then
        for name in OpenMEEG OpenMEEGMaths; do
            install_name_tool -change "${LIBGFORTRAN}" "@rpath/${GFORTRAN_NAME}" ./build/${name}/lib${name}.1.1.0.dylib
        done
    fi
    echo "::endgroup::"
    echo "::group::cmake --target package"
    cmake --build build --target package --target install --config release
    mkdir -p installers
    cp -av build/OpenMEEG-*-*.* installers/
    echo "::endgroup::"
fi

# Put DLLs where they can be found
if [[ "$PLATFORM" == 'linux'* ]]; then
    mkdir -p /output
    if [[ "$KIND" == "wheel" ]]; then
        ls -al install/lib64/*.so*
        cp -av install/lib64/*.so* /usr/local/lib/
    else
        cp -av build/OpenMEEG-*-*.* /output/
    fi
elif [[ "$PLATFORM" == 'macosx-'* ]]; then
    if [[ "$KIND" == "app" ]]; then
        otool -L $ROOT/build/OpenMEEG/libOpenMEEG.1.1.0.dylib
        otool -L $ROOT/install/lib/libgfortran*.dylib
    else
        if [[ "$PLATFORM" == 'macosx-arm64' ]]; then
            # https://matthew-brett.github.io/docosx/mac_runtime_link.html
            #cp -av $ROOT/vcpkg_installed/arm64-osx-release-11.0/lib/libomp* $ROOT/install/lib/
            otool -L $ROOT/install/lib/libOpenMEEG.1.1.0.dylib
            # install_name_tool -change "@@HOMEBREW_PREFIX@@/opt/libomp/lib/libomp.dylib" "@loader_path/libomp.dylib" $ROOT/install/lib/libOpenMEEG.1.1.0.dylib
            # otool -L $ROOT/install/lib/libOpenMEEG.1.1.0.dylib
        fi
    fi
elif [[ "$PLATFORM" == 'win'* ]]; then
    if [[ "$KIND" == "wheel" ]]; then
        cp -av $OPENBLAS_LIB/libopenblas_v0.3.20-140-gbfd9c1b5-gcc_8_1_0.dll install/bin/
    else
        ./build_tools/install_dependency_walker.sh
    fi
fi

if [[ "$KIND" == "wheel" ]]; then
    # TODO: This is only necessary because SWIG does not work outside cmake yet,
    # and we want this on windows
    if [[ "$PLATFORM" == 'win'* ]]; then
        mv build build_nopython
    fi
    echo "ls -al $PWD:"
    ls -al
fi

echo "git status:"
git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"
