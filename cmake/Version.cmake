# Derive the version from git tags, the way setuptools_scm already does for the
# Python package, so that cutting a release is just pushing a tag.
#
# Sets OPENMEEG_VERSION (full PEP 440 string, e.g. "2.5.16" on a tag or
# "2.5.17.dev68" after it) and OPENMEEG_VERSION_NUMERIC (the X.Y.Z part, which is
# all project(VERSION) accepts). Included before project(), so nothing here may
# depend on the compiler or on CMAKE_MODULE_PATH.
#
# Sources, in order: an explicit -DOPENMEEG_VERSION=, `git describe`, and
# .git_archival.txt (export-subst fills it in in GitHub's source archives).

# Turn "<tag>[-<distance>-g<sha>]" into OPENMEEG_VERSION{,_NUMERIC}.
function(_openmeeg_version_from_describe describe)
    if (describe MATCHES "^(.+)-([0-9]+)-g[0-9a-f]+$")
        set(tag "${CMAKE_MATCH_1}")
        set(distance "${CMAKE_MATCH_2}")
    else()  # %(describe) omits the suffix when HEAD is exactly on a tag
        set(tag "${describe}")
        set(distance 0)
    endif()

    # Mirror setuptools_scm's default tag parsing: an optional "name-" prefix
    # (we have a "release-2.1" tag) and an optional "v" both get dropped.
    string(REGEX REPLACE "^[A-Za-z][A-Za-z0-9_]*-" "" tag "${tag}")
    string(REGEX REPLACE "^[vV]" "" tag "${tag}")
    if (NOT tag MATCHES "^([0-9]+)\\.([0-9]+)\\.?([0-9]*)$")
        message(FATAL_ERROR "Could not parse a version out of git tag '${tag}' (from '${describe}')")
    endif()
    set(major "${CMAKE_MATCH_1}")
    set(minor "${CMAKE_MATCH_2}")
    set(patch "${CMAKE_MATCH_3}")
    if (patch STREQUAL "")
        set(patch 0)
    endif()

    if (distance EQUAL 0)
        set(version "${major}.${minor}.${patch}")
    else()
        # setuptools_scm's guess-next-dev: the tag is released, so we are working
        # towards the next patch release. local_scheme is no-local-version in
        # pyproject.toml, hence no "+g<sha>" here either.
        math(EXPR patch "${patch}+1")
        set(version "${major}.${minor}.${patch}.dev${distance}")
    endif()

    set(OPENMEEG_VERSION "${version}" PARENT_SCOPE)
    set(OPENMEEG_VERSION_NUMERIC "${major}.${minor}.${patch}" PARENT_SCOPE)
endfunction()

function(_openmeeg_describe_head source_dir out_var)
    set(${out_var} "" PARENT_SCOPE)

    # .git is a file, not a directory, in worktrees and submodules.
    find_program(GIT_EXECUTABLE NAMES git)
    if (GIT_EXECUTABLE AND EXISTS "${source_dir}/.git")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags --long --match "*[0-9]*"
            WORKING_DIRECTORY "${source_dir}"
            OUTPUT_VARIABLE describe
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_VARIABLE git_error
            RESULT_VARIABLE git_result
        )
        if (git_result EQUAL 0)
            set(${out_var} "${describe}" PARENT_SCOPE)
            return()
        endif()
        message(WARNING "git describe failed, falling back to .git_archival.txt: ${git_error}")
    endif()

    if (EXISTS "${source_dir}/.git_archival.txt")
        file(STRINGS "${source_dir}/.git_archival.txt" describe_line REGEX "^describe-name:")
        string(REGEX REPLACE "^describe-name:[ \t]*" "" describe "${describe_line}")
        # Unsubstituted in a plain checkout: export-subst only runs for archives.
        if (describe AND NOT describe MATCHES "^\\$Format:")
            set(${out_var} "${describe}" PARENT_SCOPE)
        endif()
    endif()
endfunction()

set(OPENMEEG_VERSION "" CACHE STRING "OpenMEEG version; derived from git tags when empty")

if (OPENMEEG_VERSION)
    _openmeeg_version_from_describe("${OPENMEEG_VERSION}")
else()
    _openmeeg_describe_head("${CMAKE_CURRENT_LIST_DIR}/.." describe)
    if (NOT describe)
        message(FATAL_ERROR
            "Cannot determine the OpenMEEG version: no usable git checkout and no "
            "substituted .git_archival.txt. Pass -DOPENMEEG_VERSION=X.Y.Z instead.")
    endif()
    _openmeeg_version_from_describe("${describe}")
endif()

message(STATUS "OpenMEEG version ${OPENMEEG_VERSION} (project version ${OPENMEEG_VERSION_NUMERIC})")
