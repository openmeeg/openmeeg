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

IF (WIN32)
    TARGET_LINK_LIBRARIES(OpenMEEGMaths ${LAPACK_LIBRARIES})
ENDIF()

SET_TARGET_PROPERTIES(OpenMEEGMaths PROPERTIES
                      VERSION 1.1.0
                      SOVERSION 1
                      CLEAN_DIRECT_OUTPUT 1)

INSTALL(TARGETS OpenMEEGMaths
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# install headers
FILE(GLOB MATLIB_HEADERS "${CMAKE_SOURCE_DIR}/*.{h,H}")
INSTALL(FILES ${MATLIB_HEADERS}
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)
