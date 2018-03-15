# OpenMEEG
#
# Copyright (c) INRIA 2013-2017. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

function(msinttypes_project)

    # Prepare the project

    EP_Initialisation(msinttypes USE_SYSTEM OFF BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

    # Define repository where get the sources

    if (NOT DEFINED ${ep}_SOURCE_DIR)
        set(location
            SVN_REPOSITORY "http://msinttypes.googlecode.com/svn/trunk/")
    endif()

    # set compilation flags

    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
    endif()

    set(cmake_args
        ${ep_common_cache_args}
        ${ep_optional_args}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    )

    # Add external-project

    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/patches/msinttypes.patch ${source_dir}/CMakeLists.txt
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS ${cmake_args}
    )

    # Set variable to provide infos about the project

    ExternalProject_Get_Property(${ep} install_dir)
    set(${ep}_CMAKE_FLAGS ${install_dir} PARENT_SCOPE)

    # Add custom targets

    EP_AddCustomTargets(${ep})

endfunction()
