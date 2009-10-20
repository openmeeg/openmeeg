# =============
# = Maths Lib =
# =============

ADD_LIBRARY(OpenMEEGMaths SHARED
    vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp fast_sparse_matrix.cpp MathsIO.C MatlabIO.C AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)

IF (USE_MATIO)
    TARGET_LINK_LIBRARIES(OpenMEEGMaths ${MATIO_LIBRARIES})
ENDIF()

IF (WIN32)
    TARGET_LINK_LIBRARIES(OpenMEEGMaths ${LAPACK_LIBRARIES})
ENDIF()

#SET_TARGET_PROPERTIES(OpenMEEGMaths PROPERTIES
#                      VERSION 0.1.0
#                      SOVERSION 0
#                      CLEAN_DIRECT_OUTPUT 1)

INSTALL(TARGETS OpenMEEGMaths
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
