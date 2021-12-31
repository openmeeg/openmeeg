add_library(OpenMEEG SHARED
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

set_target_properties(OpenMEEG PROPERTIES VERSION 1.1.0 SOVERSION 1 CLEAN_DIRECT_OUTPUT 1)

target_compile_definitions(OpenMEEG PUBLIC HAVE_ISNORMAL_IN_NAMESPACE_STD)
target_compile_definitions(OpenMEEG PUBLIC HAVE_ISNORMAL_IN_MATH_H)

if (VTK_FOUND)
    target_compile_definitions(OpenMEEG PUBLIC USE_VTK)
endif()

# OpenMP

if (OpenMP_FOUND)
    target_compile_definitions(OpenMEEG PUBLIC ${OPENMP_DEFINITIONS})
    target_link_libraries(OpenMEEG PUBLIC OpenMP::OpenMP_CXX)
else()
    message("OpenMP library not found. Use a compiler with OpenMP support for optimized running time." )
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
          $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # XXX: this is to include the generated file OpenMEEG_Export.h
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>         # XXX: this is to include the generated file OpenMEEGConfigure.h
          ${VTK_INCLUDE_DIRS}
          ${GIFTI_INCLUDE_DIRS}
      PRIVATE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEG/src>
)

target_link_libraries(OpenMEEG
      PUBLIC # XXX : this should become PRIVATE ! and transitive
          OpenMEEG::OpenMEEGMaths
          ${VTK_LIBRARIES}
          ${GIFTI_LIBRARIES}
)

add_library(OpenMEEG::OpenMEEG ALIAS OpenMEEG)

generate_export_header(OpenMEEG
    EXPORT_FILE_NAME OpenMEEG_Export.h
    STATIC_DEFINE OpenMEEG_BUILT_AS_STATIC)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS OpenMEEG EXPORT OpenMEEGConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)
