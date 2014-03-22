# OpenMEEG
#
# Copyright (c) INRIA 2013-2014. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

function(hdf5_project)

    set(ep hdf5)

    # List the dependencies of the project

    set(${ep}_dependencies ${MSINTTYPES} zlib)
      
    # Prepare the project

    EP_Initialisation(${ep} 
        USE_SYSTEM OFF 
        BUILD_SHARED_LIBS ON
    )

    if (NOT USE_SYSTEM_${ep})

        # Set directories

        EP_SetDirectories(${ep} EP_DIRECTORIES ep_dirs)

        # Define repository where get the sources

        if (NOT DEFINED ${ep}_SOURCE_DIR)
            set(location
                URL "http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.bz2"
                URL_MD5 "03ad766d225f5e872eb3e5ce95524a08")
        endif()

        # Add specific cmake arguments for configuration step of the project

        set(ep_optional_args)

        # set compilation flags

        if (UNIX)
            set(${ep}_c_flags "${${ep}_c_flags} -w")
        endif()

        set(cmake_args
            ${ep_common_cache_args}
            ${ep_optional_args}
            -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
            -DZLIB_ROOT:STRING=${zlib_DIR}
            -DCMAKE_C_FLAGS:STRING=${${ep}_c_flags}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_SHARED_LINKER_FLAGS:STRING=${${ep}_shared_linker_flags}  
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS_${ep}}
            -DBUILD_TESTING:BOOL=OFF
        )
        message(":::${cmake_args}:::")

        # Check if patch has to be applied

        ep_GeneratePatchCommand(${ep} PATCH_COMMAND)

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
        set(${ep}_DIR ${install_dir} PARENT_SCOPE)

        # Add custom targets

        EP_AddCustomTargets(${ep})

    endif()

endfunction()
