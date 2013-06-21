IF (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    INCLUDE(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
ENDIF()

SET(OPENMEEG_HEADERS
    analytics.h assemble.h conditioning.h cpuChrono.h danielsson.h DLLDefinesOpenMEEG.h domain.h forward.h gain.h geometry.h integrator.h
    interface.h inversers.h mesh.h om_utils.h operators.h options.h PropertiesSpecialized.h reader.h sensors.h triangle.h
    Triangle_triangle_intersection.h vect3.h vertex.h 
#   These files are imported from another repository.
#   Please do not update them in this repository.
    DataTag.H FileExceptions.H MeshDescription/Exceptions.H Properties.H
)

ADD_LIBRARY(OpenMEEG ${LIB_TYPE}
    assembleFerguson.cpp assembleHeadMat.cpp assembleSourceMat.cpp assembleSensors.cpp
    danielsson.cpp geometry.cpp operators.cpp sensors.cpp inversers.cpp cpuChrono.cpp ${OPENMEEG_HEADERS}
)
TARGET_LINK_LIBRARIES(OpenMEEG OpenMEEGMaths ${OPENMEEG_LIBRARIES} ${LAPACK_LIBRARIES})

SET_TARGET_PROPERTIES(OpenMEEG PROPERTIES
                         VERSION 1.1.0
                         SOVERSION 1
                         CLEAN_DIRECT_OUTPUT 1)

MACRO(OPENMEEG_EXECUTABLE EXECUTABLE)
    NEW_EXECUTABLE(${EXECUTABLE} ${ARGN} LIBRARIES OpenMEEG)
ENDMACRO()

OPENMEEG_EXECUTABLE(om_presmooth presmooth.cpp)
OPENMEEG_EXECUTABLE(om_minverser minverser.cpp)
OPENMEEG_EXECUTABLE(om_inverse inverse.cpp)
OPENMEEG_EXECUTABLE(om_gain gain.cpp)
OPENMEEG_EXECUTABLE(om_forward forward.cpp)
OPENMEEG_EXECUTABLE(om_assemble assemble.cpp)

# ================
# = INSTALLATION =
# ================

INSTALL(TARGETS OpenMEEG om_assemble om_forward om_gain om_inverse om_minverser om_presmooth
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# install headers
INSTALL(FILES ${OPENMEEG_HEADERS}
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)

# ===========
# = TESTING =
# ===========

IF(BUILD_TESTING)
    INCLUDE(TestHead.txt)

    TESTHEAD(1)
    TESTHEAD(2)

    OPTION(TEST_HEAD3 "Run tests on Head 3" OFF)

    IF (TEST_HEAD3)
        TESTHEAD(3)
    ENDIF()
ENDIF()
