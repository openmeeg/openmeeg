ADD_EXECUTABLE(presmooth presmooth.cpp)
TARGET_LINK_LIBRARIES (presmooth ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(minverser minverser.cpp)
TARGET_LINK_LIBRARIES (minverser ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(inverse inverse.cpp)
TARGET_LINK_LIBRARIES (inverse ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(gain gain.cpp)
TARGET_LINK_LIBRARIES (gain ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(forward forward.cpp)
TARGET_LINK_LIBRARIES (forward ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(assemble
                    assembleFerguson.cpp assembleMat.cpp assemble_RHSmatrix.cpp assemble_squids_EEG.cpp 
                    danielsson.cpp geometry.cpp assemble.cpp mesh3.cpp operateurs.cpp)
TARGET_LINK_LIBRARIES (assemble ${OPENMEEG_OTHER_LIBRARIES})

IF(ENABLE_PACKAGING)
    IF(UNIX)
        INSTALL(TARGETS assemble forward gain inverse minverser presmooth DESTINATION local/bin)
    ELSE(UNIX)
        INSTALL(TARGETS assemble forward gain inverse minverser presmooth DESTINATION bin)
    ENDIF(UNIX)
ENDIF(ENABLE_PACKAGING)
