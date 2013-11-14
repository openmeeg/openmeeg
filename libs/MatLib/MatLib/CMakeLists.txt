# =============
# = Maths Lib =
# =============

SET(MATLABIO "")
IF (USE_MATIO)
    SET(MATLABIO MatlabIO.C)
ENDIF()

ADD_LIBRARY(OpenMEEGMaths SHARED
    vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp fast_sparse_matrix.cpp
    MathsIO.C ${MATLABIO} AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)

IF (USE_MATIO)
    TARGET_LINK_LIBRARIES(OpenMEEGMaths ${MATIO_LIBRARIES})
ENDIF()

TARGET_LINK_LIBRARIES(OpenMEEGMaths ${LAPACK_LIBRARIES})

SET_TARGET_PROPERTIES(OpenMEEGMaths PROPERTIES
                      VERSION 1.1.0
                      SOVERSION 1
                      CLEAN_DIRECT_OUTPUT 1)

INSTALL(TARGETS OpenMEEGMaths
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# install headers
SET(MATLIB_HEADERS 
    diagmatrix.h DLLDefinesOpenMEEGMaths.h fast_sparse_matrix.h linop.h MatLibConfig.h 
    matrix.h RC.H matvectOps.h symmatrix.h sparse_matrix.h vector.h
    #   These files are imported from another repository.
    #   Please do not update them in this repository.
    AsciiIO.H BrainVisaTextureIO.H Exceptions.H IOUtils.H MathsIO.H MatlabIO.H RC.H 
    TrivialBinIO.H)
INSTALL(FILES ${MATLIB_HEADERS}
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)
