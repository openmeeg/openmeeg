# OpenMEEG
#
# Copyright (c) INRIA 2013-2014. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

function(OpenMEEG_project)
      
    # Prepare the project and list dependencies

    EP_Initialisation(OpenMEEG BUILD_SHARED_LIBS ON)
    EP_SetDependencies(${ep}_dependencies flens clapack matio ${MSINTTYPES})

    # Define repository where get the sources

    if (NOT DEFINED ${ep}_SOURCE_DIR)
        set(location GIT_REPOSITORY "${GIT_PREFIX}github.com/openmeeg/openmeeg.git")
    endif()

    # Set compilation flags

    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
    endif()

    set(cmake_args
        ${ep_common_cache_args}
        ${ep_optional_args}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_C_FLAGS:STRING=${${ep}_c_flags}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${${ep}_shared_linker_flags}  
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS_${ep}}
        -DBLASLAPACK_IMPLEMENTATION:STRING=${BLASLAPACK_IMPLEMENTATION}
        -DUSE_OMP:BOOL=${USE_OMP}
        -DENABLE_PACKAGING:BOOL=${ENABLE_PACKAGING}
        -DENABLE_PYTHON:BOOL=${ENABLE_PYTHON}
        -DBUILD_TESTING:BOOL=${BUILD_TESTING}
        -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
        ${clapack_CMAKE_FLAGS}
        ${zlib_CMAKE_FLAGS}
        ${hdf5_CMAKE_FLAGS}
        ${matio_CMAKE_FLAGS}
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(${ep} PATCH_COMMAND)

    # Add external-project

    set(tag UniLapack)
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

    ExternalProject_Get_Property(${ep} binary_dir)
    set(${ep}_DIR ${binary_dir} PARENT_SCOPE)

    # Add custom targets

    EP_AddCustomTargets(${ep} TESTS)

endfunction()
