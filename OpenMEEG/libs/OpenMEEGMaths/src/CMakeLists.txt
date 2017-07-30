# =====================
# = OpenMEEGMaths Lib =
# =====================

function(create_library libname)
    add_library(${libname} SHARED ${ARGN})
    target_link_libraries(${libname} ${LAPACK_LIBRARIES} ${matio_LIBRARIES})
    if (LAPACK_DEFINITIONS) # mainly for MKL
        target_compile_definitions(${libname} PUBLIC ${LAPACK_DEFINITIONS})
    endif()

    set_target_properties(${libname} PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          #INTERFACE_LINK_LIBRARIES "${NeededLibraries}"
                          CLEAN_DIRECT_OUTPUT 1)

    generate_export_header(${libname}
        EXPORT_FILE_NAME ${libname}_Export.h
        STATIC_DEFINE ${libname}_BUILT_AS_STATIC)

    install(TARGETS ${libname}
            EXPORT  OpenMEEGTargets
            RUNTIME DESTINATION ${INSTALL_BIN_DIR}
            ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
            LIBRARY DESTINATION ${INSTALL_LIB_DIR})

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${libname}_Export.h
            DESTINATION ${OpenMEEG_HEADER_INSTALLDIR} COMPONENT Development)
endfunction()

set(OpenMEEGMaths_SOURCES
    utils.cpp vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp fast_sparse_matrix.cpp
    MathsIO.C MatlabIO.C AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)
    
create_library(OpenMEEGMaths ${OpenMEEGMaths_SOURCES})
