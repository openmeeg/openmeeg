# OpenMEEG
#
# Copyright (c) INRIA 2013-2017. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

macro(hdf5_find_package)
#   Do nothing let OpenMEEG do the work.
endmacro()

function(hdf5_project)

    # Prepare the project and list dependencies

    EP_Initialisation(hdf5 BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    EP_SetDependencies(${ep}_dependencies ${MSINTTYPES} zlib)

    # Define repository where get the sources

    if (NOT DEFINED ${ep}_SOURCE_DIR)
        set(location GIT_REPOSITORY "${GIT_PREFIX}github.com/openmeeg/hdf5-matio.git")
        #set(location
        #    URL "http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.bz2"
        #    URL_MD5 "03ad766d225f5e872eb3e5ce95524a08")
    endif()

    # set compilation flags

    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
    endif()

    set(cmake_args
        ${ep_common_cache_args}
        ${ep_optional_args}
        -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
        ${zlib_CMAKE_FLAGS}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_C_FLAGS:STRING=${${ep}_c_flags}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${${ep}_shared_linker_flags}  
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS_${ep}}
        -DHDF5_INSTALL_LIB_DIR:STRING=${INSTALL_LIB_DIR}
        -DBUILD_TESTING:BOOL=OFF
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(${ep} PATCH_COMMAND hdf5-config.patch)

    # Add external-project

    #set(tag 1.8.14)
    set(tag master)
    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        ${PATCH_COMMAND}
        GIT_TAG ${tag}
        UPDATE_COMMAND ${GIT_BIN} pull
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS ${cmake_args}
        DEPENDS ${${ep}_dependencies}
    )

    # Set variable to provide infos about the project

    ExternalProject_Get_Property(${ep} install_dir)
    if (NOT WIN32)
        set(HDF5_CMAKE_INSTALL_DIR share/cmake)
    endif()
    set(${ep}_CMAKE_FLAGS -DHDF5_DIR:FILEPATH=${install_dir}/${HDF5_CMAKE_INSTALL_DIR} PARENT_SCOPE)

    # Add custom targets

    EP_AddCustomTargets(${ep})

endfunction()
