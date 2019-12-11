# OpenMEEGMath

add_library(OpenMEEGMaths SHARED
  src/vector.cpp src/matrix.cpp src/symmatrix.cpp src/sparse_matrix.cpp
  src/fast_sparse_matrix.cpp src/MathsIO.C src/MatlabIO.C src/AsciiIO.C
  src/BrainVisaTextureIO.C src/TrivialBinIO.C
)

set_target_properties(OpenMEEGMaths PROPERTIES
                      VERSION 1.1.0
                      SOVERSION 1
                      CLEAN_DIRECT_OUTPUT 1)

target_include_directories(OpenMEEGMaths
      PUBLIC
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEGMaths>
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
          ${Boost_INCLUDE_Dir}
          ${BLA_INCLUDE_DIR} # XXX : this should not be here, transitive properties should work
          $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # XXX: this is to include the generated file OpenMEEGMaths_Export.h
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>         # XXX: this is to include the generated file OpenMEEGConfigure.h
      PRIVATE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEGMaths/src>
          ${matio_INCLUDE_DIRS} # XXX : this should become transitive prop.
)

target_link_libraries(OpenMEEGMaths
      PUBLIC # XXX : this should become PRIVATE ! and hdf5 should come from matio
          ${LAPACK_LIBRARIES}
          HDF5::HDF5
          MATIO::MATIO
)
target_compile_definitions(OpenMEEGMaths PUBLIC ${BLA_DEFINITIONS} ${SHARED_PTR_DEFINITIONS})
add_library(OpenMEEG::OpenMEEGMaths ALIAS OpenMEEGMaths)

generate_export_header(OpenMEEGMaths
    EXPORT_FILE_NAME OpenMEEGMaths_Export.h
    STATIC_DEFINE OpenMEEGMaths_BUILT_AS_STATIC)

# create_library(OpenMEEGMaths ${OpenMEEGMaths_SOURCES})

# 'make install' to the correct locations (provided by GNUInstallDirs).

install(TARGETS OpenMEEGMaths EXPORT OpenMEEGConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)

# This are executables using OpenMEEGMath

add_subdirectory(tests)
