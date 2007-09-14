SET(OPENMEEG_HEADERS
    analytiques.h
    mainHeader.h
    cpuChrono.h mesh3.h
    danielsson.h
    om_utils.h
    MeshDataSmoother.h
    fcontainer.h
    operateurs.h
    sensors.h
    MeshReader.h
    geometry.h
    triangle.h
    integrateur.h
    vect3.h
    mesh3.h
    PropertiesSpecialized.h
    inverse.h
#   These files are imported from another repository.
#   Please do not update them in this repository.
    DataTag.H
    FileExceptions.H
    IOUtils.H
    MeshDescription/Reader.H
    MeshDescription/Normals.H
    MeshDescription/Domains.H
    MeshDescription/Exceptions.H
    Properties.H
)

ADD_LIBRARY(OpenMEEG ${LIB_TYPE}
    assembleFerguson.cpp assembleMat.cpp assemble_RHSmatrix.cpp assemble_squids_EEG.cpp
    danielsson.cpp geometry.cpp mesh3.cpp operateurs.cpp ${OPENMEEG_HEADERS} 
)
TARGET_LINK_LIBRARIES (OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(presmooth presmooth.cpp ${OPENMEEG_HEADERS})
TARGET_LINK_LIBRARIES (presmooth ${OPENMEEG_OTHER_LIBRARIES} OpenMEEG)

ADD_EXECUTABLE(minverser minverser.cpp ${OPENMEEG_HEADERS})
TARGET_LINK_LIBRARIES (minverser ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(inverse inverse.cpp ${OPENMEEG_HEADERS})
TARGET_LINK_LIBRARIES (inverse ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(gain gain.cpp ${OPENMEEG_HEADERS})
TARGET_LINK_LIBRARIES (gain ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(forward forward.cpp ${OPENMEEG_HEADERS})
TARGET_LINK_LIBRARIES (forward ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(assemble assemble.cpp)
TARGET_LINK_LIBRARIES (assemble OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

# =============
# = PACKAGING =
# =============

IF(ENABLE_PACKAGING)
    INSTALL(TARGETS OpenMEEG assemble forward gain inverse minverser presmooth
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
ENDIF(ENABLE_PACKAGING)

# ===========
# = TESTING =
# ===========

IF(BUILD_TEST)

INCLUDE(TestHead1.txt)
INCLUDE(TestHead2.txt)

OPTION(TEST_HEAD3 "Run tests on Head 3" OFF)

IF ( TEST_HEAD3 )
    INCLUDE(TestHead3.txt)
ENDIF ( TEST_HEAD3 )

ENDIF(BUILD_TEST)
