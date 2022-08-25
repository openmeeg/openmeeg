# OpenMEEGMath

include_directories(${BLA_INCLUDE_DIR})

set(OPENMEEGMATHS_SOURCES
    src/vector.cpp src/matrix.cpp src/symmatrix.cpp src/sparse_matrix.cpp
    src/fast_sparse_matrix.cpp src/MathsIO.C src/MatlabIO.C src/AsciiIO.C
    src/BrainVisaTextureIO.C src/TrivialBinIO.C)

add_compile_options(${WERROR_COMPILE_OPTION})

if (BUILD_SHARED_LIBS OR (UNIX AND NOT APPLE))
    add_library(OpenMEEGMaths SHARED ${OPENMEEGMATHS_SOURCES})
else()
    add_library(OpenMEEGMaths STATIC ${OPENMEEGMATHS_SOURCES})
endif()

set_target_properties(OpenMEEGMaths PROPERTIES
                      VERSION 1.1.0
                      SOVERSION 1
                      CLEAN_DIRECT_OUTPUT 1)

target_include_directories(OpenMEEGMaths
      PUBLIC
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEGMaths>
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/exports>
          ${BLA_INCLUDE_DIR} # XXX : this should not be here, transitive properties should work
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>         # XXX: this is to include the generated file OpenMEEGConfigure.h
      PRIVATE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEGMaths/src>
          ${matio_INCLUDE_DIRS} # XXX : this should become transitive prop.
)

target_link_libraries(OpenMEEGMaths
      PUBLIC # XXX : this should become PRIVATE ! and hdf5 should come from matio
          LAPACK::LAPACK
          BLAS::BLAS
          MATIO::MATIO
          HDF5::HDF5
)
add_library(OpenMEEG::OpenMEEGMaths ALIAS OpenMEEGMaths)

generate_export_header(OpenMEEGMaths
    EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/OpenMEEGMaths_Export.h
    STATIC_DEFINE OpenMEEGMaths_BUILT_AS_STATIC)

# create_library(OpenMEEGMaths ${OpenMEEGMaths_SOURCES})

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS OpenMEEGMaths EXPORT OpenMEEGConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)
install(DIRECTORY ${CMAKE_BINARY_DIR}/exports/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenMEEG)

# This are executables using OpenMEEGMath

add_subdirectory(tests)

if(WIN32)
    # https://gitlab.kitware.com/cmake/cmake/-/issues/22406
    # install(IMPORTED_RUNTIME_ARTIFACTS BLAS::BLAS DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
    if (BUILD_SHARED_LIBS)
        string (REPLACE ";" " " BLAS_LIBRARIES_SPLIT "${BLAS_LIBRARIES}")
        message("-- Configuring to copy BLAS: " ${BLAS_LIBRARIES_SPLIT})
        add_custom_command(
            TARGET OpenMEEGMaths POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${BLAS_LIBRARIES_SPLIT}
            $<TARGET_FILE_DIR:OpenMEEGMaths>)
    endif()
endif()
