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

macro(ep_Initialisation project BUILD_SHARED_LIBS build_shared_libs_def)

    set(ep ${project})

    # Option: do we want a static or a dynamic build ?

    option(BUILD_SHARED_LIBS_${ep} "Build shared libs for ${ep}" ${build_shared_libs_def})
    mark_as_advanced(BUILD_SHARED_LIBS_${ep})

    set(${ep}_c_flags ${ep_common_c_flags})
    set(${ep}_cxx_flags ${ep_common_cxx_flags})
    set(${ep}_shared_linker_flags ${ep_common_shared_linker_flags})

    # Add PIC flag if Static build on UNIX with amd64 arch

    if (UNIX)
        if (NOT BUILD_SHARED_LIBS_${ep} AND "${CMAKE_SYSTEM_PROCESSOR}" MATCHES 64)
            set(${ep}_c_flags "${${ep}_c_flags} -fPIC")
            set(${ep}_cxx_flags "${${ep}_cxx_flags} -fPIC")
        endif()
    endif()

    if (APPLE)
        if (BUILD_SHARED_LIBS_${ep})
            set(${ep}_shared_linker_flags "${${ep}_shared_linker_flags} -headerpad_max_install_names")
        endif()
    endif()

    # Remove dependencies with other external-project if a system version is used.

    foreach(dependence ${${ep}_dependencies})
        if (USE_SYSTEM_${dependence})
            list(REMOVE_ITEM ${ep}_dependencies ${dependence})
        endif()
    endforeach()

    # Add dependencies between the target of this project and the global target from the superproject

    foreach (target ${global_targets})
        add_dependencies(${target} ${ep}-${target})
    endforeach()

    # Define a directory for each target of the project

    set(DIR_VAR_NAMES DOWNLOAD BINARY STAMP TMP INSTALL)
    set(DIR_NAMES     ""       build  stamp tmp install)

    list(LENGTH DIR_VAR_NAMES dirnum)
    math(EXPR dirnum ${dirnum}-1)
    set(dirs PREFIX ${ep})
    foreach(i RANGE ${dirnum})
        list(GET DIR_VAR_NAMES ${i} var)
        list(GET DIR_NAMES     ${i} dir)
        set(dirs ${dirs} ${var}_DIR ${ep}/${dir})
        string(TOLOWER "${var}" varname)
        set(${varname}_dir ${ep}/${dir})
    endforeach()

    # Look for and define the source directory of the project

    if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ep}/CMakeLists.txt
        OR EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ep}/configure)
        set(${ep}_SOURCE_DIR SOURCE_DIR ${CMAKE_SOURCE_DIR}/${ep})
    endif()

    set(source_dir ${CMAKE_SOURCE_DIR}/${ep})
    set(ep_dirs ${dirs} SOURCE_DIR ${source_dir})
endmacro()
