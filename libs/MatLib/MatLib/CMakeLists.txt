# =============
# = Maths Lib =
# =============

SET(OpenMEEGMaths_SRC vector.cpp matrix.cpp symmatrix.cpp sparse_matrix.cpp MathsIO.C MatlabIO.C AsciiIO.C BrainVisaTextureIO.C TrivialBinIO.C)

ADD_LIBRARY(OpenMEEGMaths        SHARED ${OpenMEEGMaths_SRC})
ADD_LIBRARY(OpenMEEGMaths-static STATIC ${OpenMEEGMaths_SRC})

TARGET_LINK_LIBRARIES(OpenMEEGMaths        ${OPENMEEG_OTHER_LIBRARIES})
TARGET_LINK_LIBRARIES(OpenMEEGMaths-static ${OPENMEEG_OTHER_LIBRARIES})

SET_TARGET_PROPERTIES(OpenMEEGMaths PROPERTIES
                      VERSION 0.1.0
                      SOVERSION 0
                      CLEAN_DIRECT_OUTPUT 1)
SET_TARGET_PROPERTIES(OpenMEEGMaths-static PROPERTIES
                      VERSION 0.1.0
                      OUTPUT_NAME OpenMEEGMaths
                      CLEAN_DIRECT_OUTPUT 1)

INSTALL(TARGETS OpenMEEGMaths OpenMEEGMaths-static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
