##############################################################################
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

# Add common variables for all external-projects

set(ep_common_c_flags 
    "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_INIT} ${ADDITIONAL_C_FLAGS}")

set(ep_common_cxx_flags 
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_INIT} ${ADDITIONAL_CXX_FLAGS}")

set(ep_common_shared_linker_flags
    "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_INIT} ${ADDITIONAL_SHARED_LINKER_FLAGS}")

set(ep_common_cache_args
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
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
include(EP_SetDependencies)
include(EP_AddCustomTargets)
include(EP_GeneratePatchCommand)

# Include specific module of each project

file(GLOB projects_modules RELATIVE ${CMAKE_SOURCE_DIR} "cmake/externals/projects_modules/*.cmake")
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

macro(subprojects)
    
    set(SAVED_CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
    set(CMAKE_INSTALL_PREFIX "" )

    foreach (project ${ARGN})

        set(use_system_def OFF)
        if (USE_SYSTEM_${project})
            set(use_system_def ON)
        endif()

        # Option: do we want use the system version ?

        option(USE_SYSTEM_${project} "Use system installed version of ${project}" ${use_system_def})

        if (USE_SYSTEM_${project})
            if (COMMAND ${project}_find_package)
                call(${project}_find_package)
            else()
                find_package(${project} REQUIRED)
            endif()
        else()
            call(${project}_project)
            if (update-${project})
                set(update_dependencies ${update_dependencies} update-${project})
            endif()
            if (build-${project})
                set(build_dependencies ${build_dependencies} build-${project} install-${project})
            endif()
            set(SUBPROJECTS ${SUBPROJECTS} "${project}")
        endif()

    endforeach()

    set(CMAKE_INSTALL_PREFIX ${SAVED_CMAKE_INSTALL_PREFIX})

    #   Add global targets.

    add_custom_target(update
        DEPENDS ${update_dependencies}
        COMMAND echo All project have been updated. 
        && echo To build them, run the build target: 'cmake --build . --target build'
    )

    add_custom_target(build DEPENDS ${build_dependencies})

    if (BUILD_TESTING)
        add_custom_target(check DEPENDS ${TEST_TARGETS})
    endif()

    foreach (i ${INSTALL_DIRS})
        install(DIRECTORY ${i}/ USE_SOURCE_PERMISSIONS DESTINATION ${CMAKE_INSTALL_PREFIX})
    endforeach()

    add_custom_target(clean DEPENDS ${CLEAN_TARGETS})
endmacro()
