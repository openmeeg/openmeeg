# =============
# = Maths Lib =
# =============

SET(OpenMEEGMaths_SRC vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp fast_sparse_matrix.cpp MathsIO.C MatlabIO.C AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)

ADD_LIBRARY(OpenMEEGMaths ${LIB_TYPE} ${OpenMEEGMaths_SRC})

IF ( APPLE AND USE_ATLAS )
    TARGET_LINK_LIBRARIES(OpenMEEGMaths "-framework vecLib")
ELSE ( APPLE AND USE_ATLAS )
    TARGET_LINK_LIBRARIES(OpenMEEGMaths ${OPENMEEG_OTHER_LIBRARIES})
ENDIF ( APPLE AND USE_ATLAS )

SET_TARGET_PROPERTIES(OpenMEEGMaths PROPERTIES
                      VERSION 0.1.0
                      SOVERSION 0
                      CLEAN_DIRECT_OUTPUT 1)

INSTALL(TARGETS OpenMEEGMaths
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
