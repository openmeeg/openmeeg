# Overlay port: identical to upstream vcpkg ports/libaec at tag 2026.06.24,
# except that the source is fetched from GitHub instead of gitlab.dkrz.de.
#
# Rationale: gitlab.dkrz.de returns HTTP 429 (rate limited) for roughly half of
# all requests, which intermittently breaks CI. Upstream vcpkg has so far only
# corrected the GitLab project path (microsoft/vcpkg#51872, #52793), which does
# not address the availability problem. See the discussion in
# https://github.com/microsoft/vcpkg/issues/51921 -- there is not yet an
# upstream PR switching the port to vcpkg_from_github, so this overlay is still
# needed. Drop it once one lands.
#
# Deutsches-Klimarechenzentrum/libaec is DKRZ's own GitHub repository, listed by
# the GitLab README as an official source location, and is what HDF5 itself
# switched to in 2023 (HDFGroup/hdf5@2ca2a30). It is used in preference to the
# legacy https://github.com/MathisRosenhauer/libaec URL, which only reaches it
# via a GitHub repo-transfer redirect; both serve a byte-identical v1.1.6
# tarball, so the SHA512 below is the same either way.
#
# Everything below the download block is verbatim upstream, so refreshing this
# port against a newer vcpkg tag should be a small, obvious diff.
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO Deutsches-Klimarechenzentrum/libaec
    REF "v${VERSION}"
    SHA512 76df7501d1b7d91a43b525ba828f092f18d83f8ab09a9331e5758f93942a9758ad580baca8f9316b92a98639bde2e23cacbc2f33f52d0dd98ce7efe412cf43cd
    HEAD_REF main
)

string(COMPARE EQUAL "${VCPKG_LIBRARY_LINKAGE}" "static" BUILD_STATIC)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DBUILD_STATIC_LIBS=${BUILD_STATIC}
        -Dlibaec_INSTALL_CMAKEDIR=share/${PORT}
)
vcpkg_cmake_install()
vcpkg_copy_pdbs()
vcpkg_cmake_config_fixup()
vcpkg_replace_string("${CURRENT_PACKAGES_DIR}/share/libaec/libaec-config.cmake"
    "if(libaec_USE_STATIC_LIBS)"
    "if(\"${BUILD_STATIC}\") # forced by vcpkg"
)

# Compatibility with user's CMake < 3.18 (vcpkg claims support for >= 3.16):
# Make imported targets global so that libaec-config.cmake can create ALIAS targets.
set(_target_file "libaec_shared-targets")
if(BUILD_STATIC)
    set(_target_file "libaec_static-targets")
endif()
file(READ "${CURRENT_PACKAGES_DIR}/share/libaec/${_target_file}.cmake" libaec_targets)
string(REGEX REPLACE " (SHARED|STATIC) IMPORTED" " \\1 IMPORTED \${libaec_maybe_global}" libaec_targets "${libaec_targets}")
file(WRITE "${CURRENT_PACKAGES_DIR}/share/libaec/${_target_file}.cmake" "set(libaec_maybe_global \"\")
if(CMAKE_VERSION VERSION_LESS 3.18)
    set(libaec_maybe_global \"GLOBAL\")
endif()
${libaec_targets}
"
)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

file(INSTALL "${CURRENT_PORT_DIR}/usage" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE.txt")
