#------------------------------------------------------------
# VTK library
#------------------------------------------------------------

option(USE_VTK "Use VTK" OFF)

if (USE_VTK)
    set(REQUIRED "QUIET")
    if(CMAKE_PROJECT_NAME STREQUAL "OpenMEEG")
        set(REQUIRED "REQUIRED")
    endif()

    find_package(VTK ${REQUIRED} COMPONENTS vtkIOXML vtkIOLegacy NO_MODULE PATHS ${VTK_DIR})
    if (VTK_FOUND)
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            add_compile_options(-Wno-inconsistent-missing-override)
        endif()
        set(OPENMEEG_OTHER_LIBRARY_DIRECTORIES ${VTK_LIBRARY_DIRS} ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES})
        set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${VTK_INCLUDE_DIRS})
        set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${VTK_LIBRARIES})
    else()
        message("VTK not found, we will download and build it")
        set(USE_SYSTEM_VTK False CACHE BOOL "Use the VTK from the system" FORCE)
    endif()
    unset(REQUIRED)
endif()
