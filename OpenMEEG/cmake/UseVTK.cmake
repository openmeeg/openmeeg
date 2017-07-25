#------------------------------------------------------------
# VTK library
#------------------------------------------------------------

option(USE_VTK "Use VTK" OFF)

if (USE_VTK)
    set(REQUIRED "QUIET")
    if(CMAKE_PROJECT_NAME STREQUAL "OpenMEEG")
        set(REQUIRED "REQUIRED")
    endif()

    # what components do we want:
    set(VTK_FIND_COMPONENTS vtkIOXML vtkIOLegacy)
    mark_as_advanced(VTK_FIND_COMPONENTS)

    find_package(VTK ${REQUIRED} COMPONENTS ${VTK_FIND_COMPONENTS} NO_MODULE PATHS ${VTK_DIR})
    if (VTK_FOUND)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            add_compile_options(-Wno-inconsistent-missing-override)
        endif()
        if (NOT VTK_LIBRARY_DIRS)
            # hack because else it is not defined
            set(VTK_LIBRARY_DIRS ${VTK_DIR}/../..)
        endif()
        list(APPEND OpenMEEG_OTHER_LIBRARY_DIRS ${VTK_LIBRARY_DIRS})
        list(APPEND OpenMEEG_OTHER_INCLUDE_DIRS ${VTK_INCLUDE_DIRS})
        list(APPEND OpenMEEG_DEPENDENCIES VTK)
        set(CMAKE_MSVCIDE_RUN_PATH ${VTK_RUNTIME_LIBRARY_DIRS} ${CMAKE_MSVCIDE_RUN_PATH}) # specially for windows
    else() # in case we are in the SuperProject :
        message("VTK not found, we will download and build it")
        set(USE_SYSTEM_VTK False CACHE BOOL "Use the VTK from the system" FORCE)
    endif()

    unset(REQUIRED)
endif()
