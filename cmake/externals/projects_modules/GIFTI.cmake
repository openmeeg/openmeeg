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

function(GIFTI_project)

    # Prepare the project and list dependencies

    EP_Initialisation(GIFTI BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    EP_SetDependencies(${ep}_dependencies NIFTI)

    # Define repository where get the sources
    if (NOT DEFINED ${ep}_SOURCE_DIR)
        #set(location GIT_REPOSITORY "git://GIFTI.org/GIFTI.git" GIT_TAG ${tag})
        set(location
            URL "https://www.nitrc.org/frs/download.php/2262/gifticlib-1.0.9.tgz"
            URL_MD5 "4ac6342ab316136a9964ec752e384cf2")
    endif()

    # Add specific cmake arguments for configuration step of the project

    # set compilation flags
    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
        set(${ep}_cxx_flags "${${ep}_cxx_flags} -w")
        set(unix_additional_args -DGIFTI_USE_NVCONTROL:BOOL=ON)
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
        ${NIFTI_CMAKE_FLAGS}
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(GIFTI GIFTI_PATCH_COMMAND)

    # Add external-project

    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        ${GIFTI_PATCH_COMMAND}
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
