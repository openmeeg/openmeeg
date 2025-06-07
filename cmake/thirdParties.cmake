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

        if (VTK_VERSION VERSION_GREATER_EQUAL 9.1.0)
            option(ENABLE_VTK_BUG_WORKAROUND "Add an interposer library to workaround the vtk/expat XML reading bug.")
        endif()

        set(VTK_BUG_WORKAROUND_LIB)
        if (ENABLE_VTK_BUG_WORKAROUND)
            message("Enabling VTK bug workaround")
            set(vtkXMLWorkaround_SOURCES)
            foreach(source vtkXMLDataParser.cxx vtkXMLDataParser.h vtkXMLDataHeaderPrivate.h)
                set(vtkXMLWorkaround_SOURCES ${vtkXMLWorkaround_SOURCES} ${CMAKE_SOURCE_DIR}/vtkXMLWorkaround/${source})
            endforeach()

            add_library(vtkXMLWorkaround SHARED ${vtkXMLWorkaround_SOURCES})
            target_compile_definitions(vtkXMLWorkaround PRIVATE VTK_ABI_NAMESPACE_BEGIN VTK_ABI_NAMESPACE_END)
            target_include_directories(vtkXMLWorkaround PUBLIC ${VTK_INCLUDE_DIRS})
            target_link_options(vtkXMLWorkaround BEFORE PUBLIC "LINKER:-z,interpose")
            target_link_libraries(vtkXMLWorkaround PUBLIC ${VTK_LIBRARIES})
            set(VTK_BUG_WORKAROUND_LIB vtkXMLWorkaround)
            set(VTK_LIBRARIES vtkXMLWorkaround ${VTK_LIBRARIES})
            install(TARGETS vtkXMLWorkaround
                    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
                    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
                    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
        endif()
    endif()
endif()
