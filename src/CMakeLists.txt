SET(OPENMEEG_HEADERS
    analytics.h assemble.h cpuChrono.h mesh3.h danielsson.h om_utils.h operators.h sensors.h MeshReader.h geometry.h triangle.h integrator.h vect3.h
    mesh3.h PropertiesSpecialized.h inversers.h
#   These files are imported from another repository.
#   Please do not update them in this repository.
    DataTag.H FileExceptions.H IOUtils.H MeshDescription/Reader.H MeshDescription/Normals.H MeshDescription/Domains.H MeshDescription/Exceptions.H
    Properties.H
)

ADD_LIBRARY(OpenMEEG ${LIB_TYPE}
    assembleFerguson.cpp assembleHeadMat.cpp assembleSourceMat.cpp assembleSensors.cpp
    danielsson.cpp geometry.cpp mesh3.cpp operators.cpp sensors.cpp inversers.cpp cpuChrono.cpp ${OPENMEEG_HEADERS}
)
TARGET_LINK_LIBRARIES(OpenMEEG OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}" ${LAPACK_LIBRARIES})

MACRO(OPENMEEG_EXECUTABLE EXECUTABLE)
    NEW_EXECUTABLE(${EXECUTABLE} ${ARGN} LIBRARIES OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")
ENDMACRO()

OPENMEEG_EXECUTABLE(om_presmooth presmooth.cpp ${OPENMEEG_HEADERS})
OPENMEEG_EXECUTABLE(om_minverser minverser.cpp ${OPENMEEG_HEADERS})
OPENMEEG_EXECUTABLE(om_inverse inverse.cpp ${OPENMEEG_HEADERS})
OPENMEEG_EXECUTABLE(om_gain gain.cpp ${OPENMEEG_HEADERS})
OPENMEEG_EXECUTABLE(om_forward forward.cpp ${OPENMEEG_HEADERS})
OPENMEEG_EXECUTABLE(om_assemble assemble.cpp)
OPENMEEG_EXECUTABLE(om_internalPotential internalPotential.cpp ${OPEMEEG_HEADERS})

# ================
# = INSTALLATION =
# ================

INSTALL(TARGETS OpenMEEG om_assemble om_forward om_gain om_inverse
    om_minverser om_presmooth om_internalPotential
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# ===========
# = TESTING =
# ===========

IF(BUILD_TEST)
    INCLUDE(TestHead1.txt)
    INCLUDE(TestHead2.txt)

    OPTION(TEST_HEAD3 "Run tests on Head 3" OFF)

    IF (TEST_HEAD3)
        INCLUDE(TestHead3.txt)
    ENDIF (TEST_HEAD3)
ENDIF(BUILD_TEST)
