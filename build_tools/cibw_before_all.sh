#!/bin/bash

# Build and install (locally) OpenMEEG to prepare for SWIG-building the Python
# bindings separately.
#
# Run by cibuildwheel as CIBW_BEFORE_ALL for two different jobs:
# - KIND=wheel (cibuildwheel.yml)      -- Python wheels, no apps
# - KIND=app   (cibuildwheel_apps.yml) -- C++ apps and binary installers

set -eo pipefail

if [[ "$1" == "" ]]; then
    echo "Usage: $0 <PROJECT_PATH> [wheel|app]"
    exit 1
fi
KIND=${2:-wheel}
if [[ "$KIND" != "wheel" ]] && [[ "$KIND" != "app" ]]; then
    echo "Usage: $0 <PROJECT_PATH> [wheel|app]"
    echo "Got KIND=\"$KIND\", should be \"wheel\" or \"app\""
    exit 1
fi
cd $1
ROOT=$(pwd)
if [ -z "$RUNNER_OS" ]; then
    echo "RUNNER_OS is empty, it must be present in the environment"
    exit 1
fi

# Fail if we dirtied the checkout (unless explicitly allowed)
check_porcelain() {
    git status --porcelain --untracked-files=no
    test -z "$(git status --porcelain --untracked-files=no)" || test "$CHECK_PORCELAIN" == "false"
}

PLATFORM=$(python -c "import platform; print(f'{platform.system()}-{platform.machine()}')")

# Single source of truth for per-platform settings, so that everything below is
# a simple string compare on OS rather than repeated pattern matching.
# VCPKG_TRIPLET names a file in build_tools/vcpkg_triplets/.
# PACKAGE_ARCH_SUFFIX keeps installer filenames distinct per arch, since the
# release job in cibuildwheel_apps.yml merges every platform's artifacts into
# one directory (Windows builds only one arch, so it needs no suffix).
# BLAS_BACKEND is OpenBLAS everywhere except Darwin-arm64, which uses the
# system Accelerate framework (nothing to download, nothing to vendor).
BLAS_BACKEND="OpenBLAS"
case "$PLATFORM" in
    Linux-x86_64)
        OS=linux
        VCPKG_TRIPLET="x64-linux"
        LIB_OUTPUT_DIR="$ROOT/install/lib64"
        PACKAGE_ARCH_SUFFIX="_x86_64"
        ;;
    Linux-aarch64)
        OS=linux
        VCPKG_TRIPLET="arm64-linux"
        LIB_OUTPUT_DIR="$ROOT/install/lib64"
        PACKAGE_ARCH_SUFFIX="_aarch64"
        ;;
    Darwin-x86_64)
        OS=macos
        VCPKG_TRIPLET="x64-osx-release-1015"
        LIB_OUTPUT_DIR="$ROOT/install/lib"
        MACOS_MIN_VER="10.15"
        PACKAGE_ARCH_SUFFIX="_Intel"
        ;;
    Darwin-arm64)
        OS=macos
        VCPKG_TRIPLET="arm64-osx-release-133"
        LIB_OUTPUT_DIR="$ROOT/install/lib"
        MACOS_MIN_VER="13.3"  # first release with Accelerate's modern BLAS/LAPACK
        PACKAGE_ARCH_SUFFIX="_M1"
        BLAS_BACKEND="Accelerate"
        ;;
    Windows-AMD64)
        OS=windows
        VCPKG_TRIPLET="x64-windows-release-static"
        LIB_OUTPUT_DIR="$ROOT/install/bin"
        PACKAGE_ARCH_SUFFIX=""
        ;;
    *)
        echo "Unknown PLATFORM=\"$PLATFORM\""
        exit 1
        ;;
esac
echo "Using project root \"$ROOT\" on RUNNER_OS=\"${RUNNER_OS}\" PLATFORM=\"${PLATFORM}\" to set up KIND=\"$KIND\""

if [[ "$BLAS_BACKEND" == "OpenBLAS" ]]; then
    echo "::group::scipy-openblas64"
    # >=0.3.31 is required: 0.3.30.x and earlier shipped a lapacke.h whose symbols
    # were not SYMBOLPREFIX-mangled (OpenMathLib/OpenBLAS#5493, fixed upstream),
    # which used to have to be patched up with sed here.
    #
    # scipy-openblas64 is the ILP64 (64-bit integer) build; see
    # OpenMEEGMathsOpenBLASConfig.h and USE_SCIPY_OPENBLAS for the consuming side.
    python -m pip install "scipy-openblas64>=0.3.31"
    OPENBLAS_INCLUDE=$(python -c "import scipy_openblas64; print(scipy_openblas64.get_include_dir())")
    OPENBLAS_LIB_DIR=$(python -c "import scipy_openblas64; print(scipy_openblas64.get_lib_dir())")  # somewhere like "/absolute/path/to/site-packages/scipy_openblas64/lib"
    OPENBLAS_LIB_NAME=$(python -c "import scipy_openblas64; print(scipy_openblas64.get_library())")  # typically, "libscipy_openblas64_"
    echo "OPENBLAS_INCLUDE=\"$OPENBLAS_INCLUDE\""
    echo "OPENBLAS_LIB_DIR=\"$OPENBLAS_LIB_DIR\""
    ls -al "$OPENBLAS_INCLUDE" "$OPENBLAS_LIB_DIR"
    CMAKE_PREFIX_PATH_OPT="-DCMAKE_PREFIX_PATH=$OPENBLAS_LIB_DIR/cmake/openblas"
    echo "::endgroup::"
else
    echo "Using system $BLAS_BACKEND, no BLAS/LAPACK to download"
fi

check_porcelain

echo "::group::pip"
pip install --upgrade cmake "swig>=4.2" ninja
echo "cmake version: $(cmake --version | head -n 1)"
echo "ninja version: $(ninja --version | head -n 1)"
echo "swig version:  $(swig -version | head -n 1)"
echo "::endgroup::"

# Third-party C dependencies. macOS and Windows take hdf5, matio, zlib and
# libaec from vcpkg. Linux is deliberately different: HDF5 comes prebuilt in the
# h5py manylinux image (see manylinux-*-image in pyproject.toml), so matio is
# built here against that HDF5. Routing matio through vcpkg on Linux is not an
# option, because vcpkg's matio[mat73] hard-depends on vcpkg's own hdf5 and we
# would end up building and linking a second copy of it.
if [[ "$OS" == "linux" ]]; then
    export CMAKE_GENERATOR=Ninja
    echo "::group::yum"
    rpm --import https://repo.almalinux.org/almalinux/RPM-GPG-KEY-AlmaLinux
    yum -y install epel-release
    yum -y install curl zip unzip tar wget zlib-devel
    echo "::endgroup::"
    echo "::group::matio"
    set -x
    MATIO_VERSION=1.5.30  # keep in sync with vcpkg's matio port when practical
    wget https://github.com/tbeu/matio/releases/download/v${MATIO_VERSION}/matio-${MATIO_VERSION}.tar.gz
    tar xvf matio-${MATIO_VERSION}.tar.gz
    pushd matio-${MATIO_VERSION}
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
    cmake --build build --config Release --target install
    popd
    set +x
    echo "::endgroup::"
    export CMAKE_CXX_FLAGS="-I$OPENBLAS_INCLUDE"
    export DISABLE_CCACHE=1
    SHARED_OPT="-DBUILD_SHARED_LIBS=OFF"
    if [[ "$KIND" == "app" ]]; then
        export VCPKG_DEFAULT_TRIPLET="$VCPKG_TRIPLET"
        source ./build_tools/setup_vcpkg_compilation.sh
        LIBDIR_OPT="-DCMAKE_INSTALL_LIBDIR=lib"
        SOS=($OPENBLAS_LIB_DIR/*.so*)
        OLD_IFS=$IFS
        IFS=';'
        SOS=("${SOS[*]}")
        IFS=$OLD_IFS
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=${SOS}"
        echo "LIBRARIES_INSTALL_OPT=\"$LIBRARIES_INSTALL_OPT\""
    fi
elif [[ "$OS" == "macos" ]]; then
    export VCPKG_DEFAULT_TRIPLET="$VCPKG_TRIPLET"
    export SYSTEM_VERSION_OPT="-DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOS_MIN_VER}"
    source ./build_tools/setup_vcpkg_compilation.sh
    # libomp can cause segfaults on macos... maybe from version conflicts with OpenBLAS, or from being too recent?
    export OPENMP_OPT="-DUSE_OPENMP=OFF"
else  # windows
    export VCPKG_DEFAULT_TRIPLET="$VCPKG_TRIPLET"
    source ./build_tools/setup_vcpkg_compilation.sh
    pip install delvewheel "pefile!=2024.8.26"
    export SYSTEM_VERSION_OPT="-DCMAKE_SYSTEM_VERSION=7"
    OPENBLAS_DLL="$OPENBLAS_LIB_DIR\\$OPENBLAS_LIB_NAME.dll"
    OPENBLAS_DLL="$(cygpath -u "$OPENBLAS_DLL")"
    echo "OPENBLAS_DLL=${OPENBLAS_DLL}"
    test -f "$OPENBLAS_DLL"
    if [[ "$KIND" == "app" ]]; then
        LIBRARIES_INSTALL_OPT="-DEXTRA_INSTALL_LIBRARIES=$(cygpath -m ${OPENBLAS_DLL})"
    fi
fi

export PYTHON_OPT="-DENABLE_PYTHON=OFF"
export BLA_IMPL="-DBLA_IMPLEMENTATION=$BLAS_BACKEND"  # cmake_configure.sh wants the whole flag
if [[ "$BLAS_BACKEND" == "OpenBLAS" ]]; then
    export BLAS_LIBRARIES_OPT="-DUSE_SCIPY_OPENBLAS=ON"
fi
export WERROR_OPT="-DENABLE_WERROR=ON"
if [[ "${KIND}" == "wheel" ]]; then
    APP_OPT="-DENABLE_APPS=OFF"
else
    export BLA_STATIC_OPT="-DBLA_STATIC=ON"
    APP_OPT="-DENABLE_APPS=ON"
    PACKAGE_ARCH_OPT="-DPACKAGE_ARCH_SUFFIX=$PACKAGE_ARCH_SUFFIX"
fi
./build_tools/cmake_configure.sh -DCMAKE_WARN_DEPRECATED=FALSE -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=$ROOT/install -DCMAKE_INSTALL_UCRT_LIBRARIES=TRUE ${CMAKE_PREFIX_PATH_OPT} ${SHARED_OPT} ${BLAS_LIBRARIES_OPT} ${APP_OPT} ${LIBDIR_OPT} ${LIBRARIES_INSTALL_OPT} ${PACKAGE_ARCH_OPT}
echo "::group::cmake --build"
cmake --build build --target install --target package --config release
echo "::endgroup::"

if [[ "$KIND" == "app" ]]; then
    echo "::group::Copying installers"
    mkdir -p installers
    cp -av build/OpenMEEG-*-*.* installers/
    echo "::endgroup::"
fi
if [[ "$BLAS_BACKEND" == "OpenBLAS" ]]; then
    cp -av "$OPENBLAS_INCLUDE" "$ROOT/install/include/OpenBLAS"  # for build step, need it somewhere we know where it is!
    test -f "$ROOT/install/include/OpenBLAS/cblas.h"
fi

# Put DLLs where they can be found
if [[ "$OS" == "linux" ]]; then
    mkdir -p /output
    if [[ "$KIND" == "wheel" ]]; then
        ls -al $LIB_OUTPUT_DIR/*.so*
        # For some reason, copying to LIB_OUTPUT_DIR doesn't work for vendoring
        # into the wheel, so copy it somewhere else it will be found.
        # TODO: Maybe someday we can/should just use the EXTRA_INSTALL_LIBRARIES
        # mechanic for this, too
        cp -av "$LIB_OUTPUT_DIR"/*.so* /usr/local/lib/
        cp -av "$OPENBLAS_LIB_DIR"/*.so* /usr/local/lib/
    else
        cp -av build/OpenMEEG-*-*.* /output/
    fi
elif [[ "$OS" == "macos" ]]; then
    if [[ "$KIND" == "wheel" ]]; then
        otool -L $LIB_OUTPUT_DIR/libOpenMEEG.*.dylib
        if [[ "$BLAS_BACKEND" == "OpenBLAS" ]]; then
            cp -av "$OPENBLAS_LIB_DIR"/*.dylib* $LIB_OUTPUT_DIR/
        fi
    else
        otool -L $ROOT/build/OpenMEEG/libOpenMEEG.*.dylib
    fi
else  # windows
    if [[ "$KIND" == "wheel" ]]; then
        cp -av "$OPENBLAS_DLL" $LIB_OUTPUT_DIR/
    else
        ./build_tools/install_dependency_walker.sh
    fi
fi

if [[ "$KIND" == "wheel" ]]; then
    echo "ls -al $ROOT:"
    ls -al "$ROOT"
fi

echo "git status:"
check_porcelain
