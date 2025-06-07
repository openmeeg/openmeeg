set(OPENMEEG_SOURCES
    src/assembleFerguson.cpp
    src/assembleHeadMat.cpp
    src/assembleSourceMat.cpp
    src/assembleSensors.cpp
    src/domain.cpp
    src/mesh.cpp
    src/interface.cpp
    src/danielsson.cpp
    src/geometry.cpp
    src/operators.cpp
    src/sensors.cpp
    src/mesh_ios.cpp
    src/GeometryIOs.cpp
    src/triangle.cpp
)

add_compile_options(${WERROR_COMPILE_OPTION})

if (BUILD_SHARED_LIBS OR (UNIX AND NOT APPLE))
    add_library(OpenMEEG SHARED ${OPENMEEG_SOURCES})
else()
    add_library(OpenMEEG STATIC ${OPENMEEG_SOURCES})
endif()

set_target_properties(OpenMEEG PROPERTIES VERSION 1.1.0 SOVERSION 1 CLEAN_DIRECT_OUTPUT 1)

if (VTK_FOUND)
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
        target_link_libraries(vtkXMLWorkaround PUBLIC ${VTK_LIBRARIES})
        set(VTK_BUG_WORKAROUND_LIB vtkXMLWorkaround)
        install(TARGETS vtkXMLWorkaround
                ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
                LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    target_compile_definitions(OpenMEEG PUBLIC USE_VTK)
endif()

# OpenMP

if (OpenMP_FOUND)
    target_compile_definitions(OpenMEEG PUBLIC ${OPENMP_DEFINITIONS})
    target_link_libraries(OpenMEEG PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP library not found. Use a compiler with OpenMP support for optimized running time." )
endif()

# Progress bar
option(USE_PROGRESSBAR "Use progressar to display computation progress" ON)
if (USE_PROGRESSBAR)
    target_compile_definitions(OpenMEEG PUBLIC USE_PROGRESSBAR)
endif()

target_include_directories(OpenMEEG
      PUBLIC
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG>
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/exports>
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>         # XXX: this is to include the generated file OpenMEEGConfigure.h
          ${VTK_INCLUDE_DIRS}
          ${GIFTI_INCLUDE_DIRS}
      PRIVATE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEG/src>
)

target_link_libraries(OpenMEEG
      PUBLIC # XXX : this should become PRIVATE ! and transitive
          OpenMEEG::OpenMEEGMaths
          ${VTK_BUG_WORKAROUND_LIB}
          ${VTK_LIBRARIES}
          ${GIFTI_LIBRARIES}
)

add_library(OpenMEEG::OpenMEEG ALIAS OpenMEEG)

generate_export_header(OpenMEEG
    EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/OpenMEEG_Export.h
    STATIC_DEFINE OpenMEEG_BUILT_AS_STATIC)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS OpenMEEG EXPORT OpenMEEGConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)
install(DIRECTORY ${CMAKE_BINARY_DIR}/exports/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)
install(FILES ${CMAKE_BINARY_DIR}/OpenMEEGConfigure.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)
