# =============
# = OpenMEEGMaths Lib =
# =============

function(create_library libname)
    add_library(${libname} SHARED ${ARGN})
    target_link_libraries(${libname} ${LAPACK_LIBRARIES} ${matio_LIBRARIES})

    set_target_properties(${libname} PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          INTERFACE_LINK_LIBRARIES "${NeededLibraries}"
                          CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS ${libname}
            RUNTIME DESTINATION ${INSTALL_BIN_DIR}
            ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
            LIBRARY DESTINATION ${INSTALL_LIB_DIR})
endfunction()

set(OpenMEEGMaths_SOURCES
    vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp fast_sparse_matrix.cpp
    MathsIO.C MatlabIO.C AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)
    
create_library(OpenMEEGMaths ${OpenMEEGMaths_SOURCES})
