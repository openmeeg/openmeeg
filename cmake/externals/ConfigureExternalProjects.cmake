##############################################################################
#
# medInria
#
# Copyright (c) INRIA 2013. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
################################################################################

# Add common variables for all external-projects

set(CMAKE_INSTALL_PREFIX "" )  
mark_as_advanced(CMAKE_INSTALL_PREFIX)

set(ep_common_c_flags 
    "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_INIT} ${ADDITIONAL_C_FLAGS}")

set(ep_common_cxx_flags 
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_INIT} ${ADDITIONAL_CXX_FLAGS}")

set(ep_common_shared_linker_flags
    "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_INIT} ${ADDITIONAL_SHARED_LINKER_FLAGS}")

set(ep_common_cache_args
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
)

set(gen "${CMAKE_GENERATOR}")
if (CMAKE_EXTRA_GENERATOR)
    set(gen "${CMAKE_EXTRA_GENERATOR} -G ${CMAKE_GENERATOR}")
endif()

set(GITHUB_PREFIX https://github.com/)
if(${USE_GITHUB_SSH})
    set(GITHUB_PREFIX git@github.com:)
endif()

# Include cmake modules of external-project

include(ExternalProject) 

# Include common configuration steps

include(EP_Initialisation)
include(EP_SetDirectories)
include(EP_AddCustomTargets) 
include(EP_GeneratePatchCommand) 

# Include specific module of each project

file(GLOB projects_modules RELATIVE ${CMAKE_SOURCE_DIR} 
    "cmake/externals/projects_modules/*.cmake")
foreach(module ${projects_modules})
    include(${module})
endforeach()

# Call specific module of each project

macro(call func_name)
    string(REPLACE "-" "_" func ${func_name})
    file(WRITE tmp_call.cmake "${func}()")
    include(tmp_call.cmake OPTIONAL)
    file(REMOVE tmp_call.cmake)
endmacro()

# Add custom targets update, and build to explicitly update and rebuild all.

macro(add_external_projects)
    foreach (external_project ${ARGN})
        call(${external_project}_project)
        
        if(update-${external_project})
          set(update_dependencies ${update_dependencies} update-${external_project})
        endif()
        if(build-${external_project})
          set(build_dependencies ${build_dependencies} build-${external_project})
        endif()
    endforeach()

    add_custom_target(update
        DEPENDS ${update_dependencies}
        COMMAND echo All project have been updated. 
        && echo To build them, run the build target: 'cmake --build . --target build'
    )
      
    add_custom_target(build DEPENDS ${build_dependencies})
endmacro()
