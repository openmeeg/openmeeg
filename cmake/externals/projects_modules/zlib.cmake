# OpenMEEG
#
# Copyright (c) INRIA 2013-2017. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

macro(zlib_find_package)
#   Do nothing let OpenMEEG do the work.
endmacro()

function(zlib_project)

    # Prepare the project and list dependencies

    EP_Initialisation(zlib BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    EP_SetDependencies(${ep}_dependencies "${MSINTTYPES}")
      
    # Define repository where get the sources

    if (NOT DEFINED ${ep}_SOURCE_DIR)
        set(location GIT_REPOSITORY "${GIT_PREFIX}github.com/madler/zlib.git" GIT_TAG v1.2.11)
        #set(location
        #    URL "http://zlib.net/zlib-1.2.11.tar.gz"
        #    URL_MD5 "1c9f62f0778697a09d36121ead88e08e")
    endif()

    # set compilation flags

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
        -DINSTALL_LIB_DIR:PATH=${INSTALL_LIB_DIR}
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(${ep} PATCH_COMMAND zlib.patch)

    # Add external-project

    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        ${PATCH_COMMAND}
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS ${cmake_args}
        DEPENDS ${${ep}_dependencies}
    )

    # Set variable to provide infos about the project

    ExternalProject_Get_Property(${ep} install_dir)
    set(${ep}_CMAKE_FLAGS -DZLIB_ROOT:FILEPATH=${install_dir} PARENT_SCOPE)

    # Add custom targets

    EP_AddCustomTargets(${ep})

endfunction()
