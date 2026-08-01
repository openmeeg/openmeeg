#!/bin/bash -ef

# Fetch conda-forge's llvm-openmp so it can be vendored into macOS wheels/apps,
# and point FindOpenMP at it via OpenMP_ROOT.
#
# Homebrew's libomp is built for whatever macOS the runner happens to be, so
# vendoring that would silently raise the artifact's real minimum OS above the
# tag we advertise. The conda-forge build targets macOS 11.0, so every macOS
# artifact has to advertise at least that -- see MACOS_MIN_VER in
# cibw_before_all.sh and MACOSX_DEPLOYMENT_TARGET in wrapping/python/pyproject.toml.

# Go to the repo root
DIR=$(dirname "${BASH_SOURCE[0]}")
cd $DIR/..

echo "::group::setup_libomp_macos"

LIBOMP_VERSION="22.1.8"
case "$(uname -m)" in
    arm64)
        LIBOMP_SUBDIR="osx-arm64"
        LIBOMP_BUILD="hc7d1edf_0"
        ;;
    x86_64)
        LIBOMP_SUBDIR="osx-64"
        LIBOMP_BUILD="h0d3cbff_0"
        ;;
    *)
        echo "Unknown machine \"$(uname -m)\""
        exit 1
        ;;
esac
LIBOMP_PKG="llvm-openmp-${LIBOMP_VERSION}-${LIBOMP_BUILD}.conda"
LIBOMP_URL="https://api.anaconda.org/download/conda-forge/llvm-openmp/${LIBOMP_VERSION}/${LIBOMP_SUBDIR}/${LIBOMP_PKG}"

LIBOMP_ROOT="$PWD/libomp"
rm -rf "$LIBOMP_ROOT"
mkdir -p "$LIBOMP_ROOT/download"
echo "Downloading $LIBOMP_URL"
curl -sSL --retry 3 -o "$LIBOMP_ROOT/download/$LIBOMP_PKG" "$LIBOMP_URL"

# A .conda is a zip of zstd tarballs; the payload one is named pkg-*.
unzip -q "$LIBOMP_ROOT/download/$LIBOMP_PKG" -d "$LIBOMP_ROOT/download"
tar -xf "$LIBOMP_ROOT/download"/pkg-*.tar.zst -C "$LIBOMP_ROOT"

test -f "$LIBOMP_ROOT/lib/libomp.dylib"
test -f "$LIBOMP_ROOT/include/omp.h"
otool -l "$LIBOMP_ROOT/lib/libomp.dylib" | grep -A3 LC_BUILD_VERSION | grep minos

export OpenMP_ROOT="$LIBOMP_ROOT"          # honored by find_package(OpenMP)
export LIBOMP_LIB_DIR="$LIBOMP_ROOT/lib"   # for vendoring into the artifact
echo "OpenMP_ROOT=$OpenMP_ROOT"

echo "::endgroup::"
