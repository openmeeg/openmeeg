# OpenMEEG
#
# Copyright (c) INRIA 2013-2017. All rights reserved.
# See LICENSE.txt for details.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

macro(matio_find_package)
#   Do nothing let OpenMEEG do the work.
endmacro()

function(matio_project)

    # Prepare the project and list dependencies

    EP_Initialisation(matio BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    EP_SetDependencies(${ep}_dependencies ${MSINTTYPES} hdf5 zlib)

    # Define repository where get the sources

    if (NOT DEFINED ${ep}_SOURCE_DIR)
        set(location GIT_REPOSITORY "${GIT_PREFIX}github.com/openmeeg/matio-openmeeg.git")
    endif()

    # set compilation flags

    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
    endif()

    if (MSINTTYPES)
        set(MSINTTYPES_CMAKE_ARG -DINTTYPES_INCLUDES:FILEPATH=${msinttypes_DIR}/include)
    endif()

    if (zlib_CMAKE_FLAGS)
        set(zlib_CMAKE_FLAGS -DUSE_SYSTEM_ZLIB:BOOL=OFF ${zlib_CMAKE_FLAGS})
    endif()

    set(cmake_args
        ${ep_common_cache_args}
        ${ep_optional_args}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_C_FLAGS:STRING=${${ep}_c_flags}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        ${zlib_CMAKE_FLAGS}
        ${hdf5_CMAKE_FLAGS}
        -DMAT73:BOOL=ON
        ${MSINTTYPES_CMAKE_ARG}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${${ep}_shared_linker_flags}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS_${ep}}
        -DMATIO_BUILD_TESTING:BOOL=${MATIO_BUILD_TESTING}
        -DMATLAB_TESTING:BOOL=${MATLAB_TESTING}
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(${ep} PATCH_COMMAND)

    # Add external-project

    set(tag master)
    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        GIT_TAG ${tag}
        UPDATE_COMMAND ${GIT_BIN} pull
        ${PATCH_COMMAND}
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS ${cmake_args}
        DEPENDS ${${ep}_dependencies}
    )

    # Set variable to provide infos about the project

    ExternalProject_Get_Property(${ep} install_dir)
    set(${ep}_CMAKE_FLAGS -D${ep}_DIR:FILEPATH=${install_dir} PARENT_SCOPE)

    # Add custom targets

    if (MATIO_BUILD_TESTING)
        set(MATIO_TESTING TESTS)
    endif()
    EP_AddCustomTargets(${ep} ${MATIO_TESTING})

endfunction()
