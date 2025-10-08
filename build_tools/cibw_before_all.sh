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

echo "::group::scipy-openblas32"
PLATFORM=$(python -c "import platform; print(f'{platform.system()}-{platform.machine()}')")
echo "PLATFORM=$PLATFORM"
# PLATFORM can be:
# Linux-x86_64
# Linux-aarch64
# Darwin-x86_64
# Darwin-arm64
# Windows-AMD64

python -m pip install scipy-openblas32
# fix a bug in the headers!
if [[ "$PLATFORM" == 'Darwin-'* ]]; then
    SED_OPT="-i ''"
else
    SED_OPT="-i"
fi
sed $SED_OPT "s/ LAPACKE_/ scipy_LAPACKE_/g" "$(python -c 'import scipy_openblas32; print(scipy_openblas32.get_include_dir())')/lapacke.h"
OPENBLAS_INCLUDE=$(python -c "import scipy_openblas32; print(scipy_openblas32.get_include_dir())")
echo "OPENBLAS_INCLUDE=\"$OPENBLAS_INCLUDE\""
ls -alR $OPENBLAS_INCLUDE
OPENBLAS_LIB_DIR=$(python -c "import scipy_openblas32; print(scipy_openblas32.get_lib_dir())")
# mkdir -p ./.openblas
# echo "./.openblas/scipy_openblas.pc:"
# echo $(python -c "import pathlib, scipy_openblas32; pathlib.Path('./.openblas/scipy_openblas.pc').write_text(scipy_openblas32.get_pkg_config())")
# export PKG_CONFIG_PATH="$PWD/.openblas"
# echo "PKG_CONFIG_PATH=\"$PKG_CONFIG_PATH\""
# cat $PKG_CONFIG_PATH/scipy_openblas.pc
CMAKE_PREFIX_PATH_OPT="-DCMAKE_PREFIX_PATH=$OPENBLAS_LIB_DIR/cmake/openblas"
# export OpenBLAS_DIR="$OPENBLAS_LIB_DIR/cmake/openblas"
echo "CMAKE_MODULE_PATH=\"$CMAKE_MODULE_PATH\""
echo "OPENBLAS_LIB_DIR=\"$OPENBLAS_LIB_DIR\""
ls -alR $OPENBLAS_LIB_DIR
echo "::endgroup::"

git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"

if [[ "$PLATFORM" == 'Linux-'* ]]; then
    echo "::group::yum"
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux
    yum -y install epel-release
    yum -y install hdf5-devel matio-devel curl zip unzip tar ninja-build
    echo "::endgroup::"
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export LINKER_OPT="-lpthread"
    export DISABLE_CCACHE=1
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
    if [[ "$KIND" == "app" ]]; then
        if [[ "$PLATFORM" == "Linux-x86_64" ]]; then
            export VCPKG_DEFAULT_TRIPLET="x64-linux"
        elif [[ "$PLATFORM" == "Linux-aarch64" ]]; then
            export VCPKG_DEFAULT_TRIPLET="arm64-linux"
        else
            echo "Unknown PLATFORM=\"$PLATFORM\""
            exit 1
        fi
        source ./build_tools/setup_vcpkg_compilation.sh
        LAPACK_LIBRARIES_OPT="-DLAPACK_LIBRARIES=/usr/local/lib/libopenblas.a"
        LIBDIR_OPT="-DCMAKE_INSTALL_LIBDIR=lib"
    fi
elif [[ "$PLATFORM" == 'Darwin-'* ]]; then
    if [[ "$PLATFORM" == "Darwin-x86_64" ]]; then
        VC_NAME="x64"
        MIN_VER="10.15"
    elif [[ "$PLATFORM" == "Darwin-arm64" ]]; then
        VC_NAME="arm64"
        MIN_VER="11.0"
    else
        echo "Unknown PLATFORM=\"$PLATFORM\""
        exit 1
    fi
    export VCPKG_DEFAULT_TRIPLET="${VC_NAME}-osx-release-${MIN_VER//./}"  # remove dots
    export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=${MIN_VER}"
    source ./build_tools/setup_vcpkg_compilation.sh
    # libomp can cause segfaults on macos... maybe from version conflicts with OpenBLAS, or from being too recent?
    export OPENMP_OPT="-DUSE_OPENMP=OFF"
    if [[ "$KIND" == "app" ]]; then
        if [[ "$PLATFORM" == "macosx-x86_64" ]]; then
            PACKAGE_ARCH_SUFFIX="_Intel"
        else
            PACKAGE_ARCH_SUFFIX="_M1"
        fi
        PACKAGE_ARCH_OPT="-DPACKAGE_ARCH_SUFFIX=$PACKAGE_ARCH_SUFFIX"
    fi
elif [[ "$PLATFORM" == "Windows-AMD64" ]]; then
    export VCPKG_DEFAULT_TRIPLET="x64-windows-release-static"
    export CMAKE_GENERATOR="Visual Studio 17 2022"
    source ./build_tools/setup_vcpkg_compilation.sh
    pip install delvewheel "pefile!=2024.8.26"
    export SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
    if [[ "$KIND" == "app" ]]; then
        OPENBLAS_DLL=$OPENBLAS_LIB_DIR\\$(python -c "import scipy_openblas32; print(scipy_openblas32.get_library())").dll
        OPENBLAS_DLL=$(cygpath -u $OPENBLAS_DLL)
        echo "OPENBLAS_DLL=${OPENBLAS_DLL}"
        test -f $OPENBLAS_DLL
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=$(cygpath -m ${OPENBLAS_DLL})"
    fi
else
    echo "Unknown platform: ${PLATFORM}"
    exit 1
fi
export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPLEMENTATION="OpenBLAS"
export BLAS_LIBRARIES_OPT="-DUSE_SCIPY_OPENBLAS=ON"
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
    echo "::endgroup::"
    echo "::group::cmake --target package"
    cmake --build build --target package --target install --config release
    mkdir -p installers
    cp -av build/OpenMEEG-*-*.* installers/
    echo "::endgroup::"
fi
cp -a "$OPENBLAS_INCLUDE" "${ROOT}/install/include/OpenBLAS"  # for build step, need it somewhere we know where it is!
test -f "${ROOT}/install/include/OpenBLAS/cblas.h"

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
        cp -av $OPENBLAS_DLL install/bin/
    else
        ./build_tools/install_dependency_walker.sh
    fi
fi

if [[ "$KIND" == "wheel" ]]; then
    echo "ls -al $PWD:"
    ls -al
fi

echo "git status:"
git status --porcelain --untracked-files=no
test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"
