#!/bin/bash

# Make a CPack Linux tarball self-contained, like auditwheel does for wheels:
# copy every library the apps and libraries need that is not part of the base
# system into lib/, and give the copies an $ORIGIN rpath so they find each
# other (the OpenMEEG binaries already look in $ORIGIN/../lib).

set -eo pipefail
if [[ "$1" == "" ]]; then
    echo "Usage: $0 <TARBALL>"
    exit 1
fi
TARBALL=$(realpath "$1")
test -f "$TARBALL"
if ! command -v patchelf > /dev/null; then
    pip install patchelf
fi

# Libraries every manylinux_2_28 system provides, so they are not bundled (as
# for the wheels), copied from the manylinux_2_28 "lib_whitelist" in auditwheel
# 6.8.2's auditwheel/policy/manylinux-policy.json; keep in sync with that. The
# dynamic loader is never bundled either.
SYSTEM_LIBS="libGL.so.1 libICE.so.6 libSM.so.6 libX11.so.6 libXext.so.6
    libXrender.so.1 libanl.so.1 libatomic.so.1 libc.so.6 libdl.so.2
    libexpat.so.1 libgcc_s.so.1 libglib-2.0.so.0 libgobject-2.0.so.0
    libgthread-2.0.so.0 libm.so.6 libmvec.so.1 libnsl.so.1 libpthread.so.0
    libresolv.so.2 librt.so.1 libstdc++.so.6 libutil.so.1 libz.so.1"

WORK=$(mktemp -d)
tar xzf "$TARBALL" -C "$WORK"
DIR=$(echo "$WORK"/OpenMEEG-*)
test -d "$DIR/lib"

# The ELF files to check: apps and (non-symlink) libraries
elf_files() {
    find "$DIR/bin" -type f
    find "$DIR/lib" -maxdepth 1 -type f -name "*.so*"
}

# Print "<soname> <path>" for each non-system dependency of the ELF files,
# failing on any that cannot be found
external_deps() {
    local f out
    for f in $(elf_files); do
        out=$(ldd "$f")
        if grep -q "not found" <<< "$out"; then
            echo "Missing libraries for $f:" >&2
            grep "not found" <<< "$out" >&2
            return 1
        fi
        awk '/=> \// {print $1, $3}' <<< "$out"
    done | sort -u | awk -v sys="$SYSTEM_LIBS" '
        BEGIN {n = split(sys, a); for (i = 1; i <= n; i++) skip[a[i]] = 1}
        !($1 in skip) && $1 !~ /^ld-linux/'

}

DEPS=$(external_deps)
while read -r NAME LIBPATH; do
    [[ "$NAME" == "" ]] && continue
    [[ "$LIBPATH" == "$DIR"/* ]] && continue
    echo "Bundling $LIBPATH"
    cp -L "$LIBPATH" "$DIR/lib/$NAME"
    chmod u+w "$DIR/lib/$NAME"
    patchelf --set-rpath '$ORIGIN' "$DIR/lib/$NAME"
done <<< "$DEPS"

# Everything that is not a system library must now resolve inside the bundle
DEPS=$(external_deps)
LEFT=$(awk -v d="$DIR/" 'index($2, d) != 1' <<< "$DEPS")
if [[ "$LEFT" != "" ]]; then
    echo "Libraries still resolved from outside the bundle:"
    echo "$LEFT"
    exit 1
fi
ls -al "$DIR/lib"

tar czf "$TARBALL" -C "$WORK" "$(basename "$DIR")"
rm -rf "$WORK"
