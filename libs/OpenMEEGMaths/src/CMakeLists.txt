# =============
# = OpenMEEGMaths Lib =
# =============

function(create_library libname)
    add_library(${libname} SHARED ${ARGN})

    set(NeededLibraries ${LAPACK_LIBRARIES})
    if (USE_MATIO)
        set(NeededLibraries ${matio_LIBRARIES} ${NeededLibraries})
    endif()
    target_link_libraries(${libname} ${NeededLibraries})

    set_target_properties(${libname} PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          INTERFACE_LINK_LIBRARIES "${NeededLibraries}"
                          CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS ${libname}
            ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
            LIBRARY DESTINATION ${INSTALL_LIB_DIR})
endfunction()

set(MATLABIO "")
if (USE_MATIO)
    set(MATLABIO MatlabIO.C)
endif()

set(OpenMEEGMaths_SOURCES
    vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp fast_sparse_matrix.cpp
    MathsIO.C ${MATLABIO} AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)
    
create_library(OpenMEEGMaths ${OpenMEEGMaths_SOURCES})
#create_library(OpenMEEGMathsInteractive ${OpenMEEGMaths_SOURCES})
#set_target_properties(OpenMEEGMathsInteractive PROPERTIES COMPILE_DEFINITIONS "IGNORE_NDEBUG=1")
