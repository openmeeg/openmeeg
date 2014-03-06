# =============
# = Maths Lib =
# =============

function(create_library libname)
    add_library(${libname} SHARED ${ARGN})

    if (USE_MATIO)
        target_link_libraries(${libname} ${MATIO_LIBRARIES})
    endif()

    target_link_libraries(${libname} ${LAPACK_LIBRARIES})

    set_target_properties(${libname} PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS ${libname}
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
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

# install headers
set(MATLIB_HEADERS 
    DLLDefinesOpenMEEGMaths.h fast_sparse_matrix.h linop.h MatLibConfig.h 
    matrix.h RC.H matvectOps.h symmatrix.h sparse_matrix.h vector.h
    #   These files are imported from another repository.
    #   Please do not update them in this repository.
    AsciiIO.H BrainVisaTextureIO.H Exceptions.H IOUtils.H MathsIO.H MatlabIO.H RC.H 
    TrivialBinIO.H)

install(FILES ${MATLIB_HEADERS}
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)
