include(CMakePackageConfigHelpers)

#   Offer the user the choice of overriding the installation directories

get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

set(LIBSUFFIX "")
if ("X${LIB64}" STREQUAL "XTRUE" AND NOT APPLE)
    set(LIBSUFFIX 64)
endif()

set(INSTALL_LIB_DIR     lib${LIBSUFFIX} CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin             CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include         CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR    share           CACHE PATH "Installation directory for data files")
set(INSTALL_CMAKE_DIR   share/${CMAKE_PROJECT_NAME}/cmake CACHE PATH "Installation directory for CMake files")

if (WIN32)
    set(INSTALL_DATA_DIR ${CMAKE_PROJECT_NAME} CACHE PATH "Installation directory for data files")
    set(INSTALL_CMAKE_DIR cmake CACHE PATH "Installation directory for CMake files")
endif()

mark_as_advanced(INSTALL_LIB_DIR INSTALL_BIN_DIR INSTALL_INCLUDE_DIR INSTALL_DATA_DIR)
 
macro(sub_directories)
    foreach (dir ${ARGN})
        add_subdirectory(${dir})
    endforeach()
endmacro()

macro(PARSE_ARGUMENTS_GENERATE LIST_VARS DEFAULT_VAR)
    foreach(var ${LIST_VARS})
        unset(${var})
    endforeach ()

    set(CURRENT_VAR ${DEFAULT_VAR})
    foreach (arg ${ARGN})
        set(skip_this_arg FALSE)
        foreach(var ${LIST_VARS})
            if (${arg} STREQUAL ${var})
                set(CURRENT_VAR ${var})
                set(skip_this_arg TRUE)
                break()
            endif()
        endforeach ()
        if (NOT skip_this_arg)
            set(${CURRENT_VAR} ${${CURRENT_VAR}} ${arg})
        endif()
    endforeach()
endmacro()

macro(Dependencies ConfigName)
    PARSE_ARGUMENTS_GENERATE("OPTIONAL" "DEFAULT" ${ARGN})
    set(DepFileName ${CMAKE_CURRENT_BINARY_DIR}/${ConfigName}Dependencies.cmake)
    set(DEPLIBS)
    set(DEPINCS)
    set(FindDep) # <- contains the list of FindDep.cmake needed to be installed
    message("Dependencies for ${ConfigName}")
    if (NOT "${DEFAULT}" STREQUAL "" OR NOT "${OPTIONAL}" STREQUAL "")
        configure_file(${PROJECT_SOURCE_DIR}/cmake/FindSoftware.cmake ${DepFileName} COPYONLY)
        file(APPEND ${DepFileName} "\n")
        file(APPEND ${DepFileName} "set(FIND_DEBUG_MODE 1)\n")
        set(REQUIRED " REQUIRED")
        foreach(arg ${DEFAULT} NOT_REQUIRED ${OPTIONAL})
            if ("${arg}" STREQUAL "NOT_REQUIRED")
                set(REQUIRED "")
            else()
                file(APPEND ${DepFileName} "find(package ${arg} COMPONENTS ${${arg}_FIND_COMPONENTS} PATHS \${${arg}_DIR} NO_DEFAULT_PATH QUIET CONFIG)\n")
                file(APPEND ${DepFileName} "find(package ${arg} COMPONENTS ${${arg}_FIND_COMPONENTS} MODULE QUIET)\n")
                file(APPEND ${DepFileName} "find(package ${arg}${REQUIRED})\n\n")
                set(DEPLIBS "${DEPLIBS} \${${arg}_LIBRARIES}")
                set(DEPINCS "${DEPINCS} \${${arg}_INCLUDE_DIRS} \${${arg}_INCLUDE_DIR}")
                if (EXISTS ${PROJECT_SOURCE_DIR}/cmake/Find${arg}.cmake)
                    list(APPEND FindDep ${PROJECT_SOURCE_DIR}/cmake/Find${arg}.cmake)
                endif()
            endif()
        endforeach()
        include(${DepFileName})
        install(FILES ${DepFileName} ${FindDep} DESTINATION ${INSTALL_CMAKE_DIR})
    endif()
endmacro()

function(GenerateConfigFile ConfigName)

    PARSE_ARGUMENTS_GENERATE("LIBRARIES;INCLUDE_DIRS" "DEFAULT" ${ARGN})

    #   Make relative paths absolute (needed later on)
    foreach(p LIB BIN INCLUDE DATA)
        set(var INSTALL_${p}_DIR)
        if(NOT IS_ABSOLUTE "${${var}}")
            set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
        endif()
    endforeach()

    #   Creating files for companion projects
    set(LIBRARIES "${LIBRARIES}")
    set(version "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

    #   Create a XXXConfig.cmake file for the use from the build tree
    foreach (dir ${INCLUDE_DIRS})
        set(I_D ${I_D} "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
    endforeach()
    set(${ConfigName}_INCLUDE_DIRS ${I_D} "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/include")

    set(${ConfigName}_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}")

    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/XXXConfig.cmake.in "${ConfigName}Config.cmake"
                                  INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}
                                  PATH_VARS INSTALL_INCLUDE_DIR 
                                            INSTALL_DATA_DIR
                                            INSTALL_LIB_DIR
                                            INSTALL_CMAKE_DIR
                                  )

    write_basic_package_version_file("${ConfigName}ConfigVersion.cmake"
                                     VERSION ${version}
                                     COMPATIBILITY SameMajorVersion)
    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/${ConfigName}Config.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/${ConfigName}ConfigVersion.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/Use${ConfigName}.cmake"
            DESTINATION ${INSTALL_CMAKE_DIR})
endfunction()
