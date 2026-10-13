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

# glibc, the base GCC runtime and zlib, which every manylinux-compatible system
# provides (cf. auditwheel's policies); everything else gets bundled
SYSTEM_RE='^(linux-vdso|linux-gate|ld-linux.*|libc|libm|libdl|librt|libpthread|libutil|libresolv|libnsl|libstdc\+\+|libgcc_s|libz)\.so'

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
    done | sort -u | { grep -Ev "$SYSTEM_RE" || true; }
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
