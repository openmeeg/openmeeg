################################################################################
#
# OpenMEEG
#
# Copyright (c) INRIA 2013-2017. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
################################################################################

function(NIFTI_project)

    # Prepare the project and list dependencies

    EP_Initialisation(NIFTI BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    EP_SetDependencies(${ep}_dependencies zlib)

    # Define repository where get the sources
    if (NOT DEFINED ${ep}_SOURCE_DIR)
        #set(location GIT_REPOSITORY "git://NIFTI.org/NIFTI.git" GIT_TAG ${tag})
        set(location
            URL "https://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz"
            URL_MD5 "425a711f8f92fb1e1f088cbc55bea53a")
    endif()

    # Add specific cmake arguments for configuration step of the project

    # set compilation flags
    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
        set(${ep}_cxx_flags "${${ep}_cxx_flags} -w")
        set(unix_additional_args -DNIFTI_USE_NVCONTROL:BOOL=ON)
    endif()

    set(cmake_args
        ${ep_common_cache_args}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_PREFIX_PATH:PATH=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_C_FLAGS:STRING=${${ep}_c_flags}
        -DCMAKE_CXX_FLAGS:STRING=${${ep}_cxx_flags}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${${ep}_shared_linker_flags}  
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS_${ep}}
        -DBUILD_TESTING:BOOL=OFF
        ${zlib_CMAKE_FLAGS}
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(NIFTI NIFTI_PATCH_COMMAND)

    # Add external-project

    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        ${NIFTI_PATCH_COMMAND}
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS ${cmake_args}
        DEPENDS ${${ep}_dependencies}
    )
      
    # Set variable to provide infos about the project

    ExternalProject_Get_Property(${ep} install_dir)
    set(${ep}_CMAKE_FLAGS -D${ep}_DIR:FILEPATH=${install_dir} PARENT_SCOPE)

    # Add custom targets

    EP_AddCustomTargets(${ep})

endfunction()
