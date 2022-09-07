include(FindBLASImplementation)

if (BLA_STATIC)
    set(MATIO_USE_STATIC_LIBRARIES TRUE) # XXX This should be an option
endif()

find_package(matio REQUIRED)

# VTK stuff

if (USE_VTK)
    # what components do we want:
    set(VTK_COMPONENTS IOXML IOLegacy)
    if (${CMAKE_VERSION} VERSION_LESS_EQUAL "3.18")
        set(VTK_COMPONENTS vtkIOXML vtkIOLegacy)
    endif()
    mark_as_advanced(VTK_COMPONENTS)

    find_package(VTK REQUIRED COMPONENTS ${VTK_COMPONENTS})
    if (VTK_FOUND)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            add_compile_options(-Wno-inconsistent-missing-override)
        endif()
        # set(CMAKE_MSVCIDE_RUN_PATH ${VTK_RUNTIME_LIBRARY_DIRS} ${CMAKE_MSVCIDE_RUN_PATH}) # specially for windows
        message(STATUS "Found VTK, including requested VTK IO support...")
    endif()
endif()
