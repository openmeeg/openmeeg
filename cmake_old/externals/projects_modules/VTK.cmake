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

function(VTK_project)

    # Prepare the project and list dependencies

    EP_Initialisation(VTK BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    EP_SetDependencies(${ep}_dependencies zlib)

    # Define repository where get the sources

    # Set GIT_TAG to latest commit of origin/release-6.3 known to work
    set(tag 9e24f51afcaebd4fbd474e8f9e620bad8997c0a3)
    if (NOT DEFINED ${ep}_SOURCE_DIR)
        #set(location GIT_REPOSITORY "git://vtk.org/VTK.git" GIT_TAG ${tag})
        set(location
            URL "http://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz"
            URL_MD5 "daee43460f4e95547f0635240ffbc9cb")
    endif()

    # Add specific cmake arguments for configuration step of the project

    # set compilation flags
    if (UNIX)
        set(${ep}_c_flags "${${ep}_c_flags} -w")
        set(${ep}_cxx_flags "${${ep}_cxx_flags} -w")
        set(unix_additional_args -DVTK_USE_NVCONTROL:BOOL=ON)
    endif()

    set(cmake_args
        ${ep_common_cache_args}
        ${zlib_CMAKE_FLAGS}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_SKIP_INSTALL_RPATH:BOOL=ON
        -DCMAKE_SKIP_RPATH:BOOL=ON
        -DCMAKE_C_FLAGS:STRING=${${ep}_c_flags}
        -DCMAKE_CXX_FLAGS:STRING=${${ep}_cxx_flags}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${${ep}_shared_linker_flags}  
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS_${ep}}
        -DVTK_WRAP_TCL:BOOL=OFF
        -DVTK_USE_SYSTEM_ZLIB:BOOL=ON
        -DVTK_Group_Rendering:BOOL=OFF
        -DVTK_Group_StandAlone:BOOL=OFF
        -DModule_vtkIOLegacy:BOOL=ON
        -DModule_vtkIOCore:BOOL=ON
        -DModule_vtkIOXML:BOOL=ON
        -DBUILD_TESTING:BOOL=OFF 
    )

    # Check if patch has to be applied

    ep_GeneratePatchCommand(VTK VTK_PATCH_COMMAND)

    # Add external-project

    ExternalProject_Add(${ep}
        ${ep_dirs}
        ${location}
        ${VTK_PATCH_COMMAND}
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
