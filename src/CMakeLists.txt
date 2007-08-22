SET(OPENMEEG_HEADERS
    DataTag.h
    analytiques.h
    mainHeader.h
    FileExceptions.h
    cpuChrono.h mesh3.h
    IOUtils.h
    danielsson.h
    om_utils.h
    MeshDataSmoother.h
    fcontainer.h
    operateurs.h
    MeshDescriptionReader.h
    genericMatrix.h
    sensors.h
    MeshDescriptionReaderSpecialized.h
    geometry.h
    triangle.h
    Properties.h
    integrateur.h
    vect3.h
    mesh3.h
    PropertiesSpecialized.h
    inverse.h
)

ADD_LIBRARY(OpenMEEG
    assembleFerguson.cpp assembleMat.cpp assemble_RHSmatrix.cpp assemble_squids_EEG.cpp
    danielsson.cpp geometry.cpp mesh3.cpp operateurs.cpp ${OPENMEEG_HEADERS}
)

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
    INSTALL(TARGETS assemble forward gain inverse minverser presmooth DESTINATION bin)
ENDIF(ENABLE_PACKAGING)

# ===========
# = TESTING =
# ===========

IF(BUILD_TEST)

############ ASSEMBLE TEST ##############
ADD_TEST(assemble-LHS ${CMAKE_CURRENT_BINARY_DIR}/assemble -LHS ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin)

ADD_TEST(minverser-lhsmatrix ${CMAKE_CURRENT_BINARY_DIR}/minverser ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin)

ADD_TEST(assemble-RHS ${CMAKE_CURRENT_BINARY_DIR}/assemble -RHS ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin)

ADD_TEST(presmooth ${CMAKE_CURRENT_BINARY_DIR}/presmooth ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin)


############ EEG TEST ##############
ADD_TEST(assemble-x2EEGresp ${CMAKE_CURRENT_BINARY_DIR}/assemble -xToEEGresponse ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.patches ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin)

ADD_TEST(gain-EEG ${CMAKE_CURRENT_BINARY_DIR}/gain -EEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin)

#------ FORWARD -----#
ADD_TEST(forward-EEG ${CMAKE_CURRENT_BINARY_DIR}/forward -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.src ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt 0.0)

#------ INVERSE -----#
ADD_TEST(inverse-EEG-HEAT ${CMAKE_CURRENT_BINARY_DIR}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGHEAT.txt 1 1e-4 HEAT 10000 0.0)

ADD_TEST(inverse-EEG-MN ${CMAKE_CURRENT_BINARY_DIR}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGMN.txt 1 1e-4 MN 10000 0.0)

ADD_TEST(inverse-EEG-TV ${CMAKE_CURRENT_BINARY_DIR}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGTV.txt 1 1e-4 TV 10000 0.0)


############ MEG TEST ##############
ADD_TEST(assemble-x2MEGresp ${CMAKE_CURRENT_BINARY_DIR}/assemble -xToMEGresponseContrib ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin)

ADD_TEST(assemble-s2MEGresp ${CMAKE_CURRENT_BINARY_DIR}/assemble -sToMEGresponseContrib ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin)

ADD_TEST(gain-MEG ${CMAKE_CURRENT_BINARY_DIR}/gain -MEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin)

#------ FORWARD -----#
ADD_TEST(forward-MEG ${CMAKE_CURRENT_BINARY_DIR}/forward -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.src ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt 0.0)

#------ INVERSE -----#
ADD_TEST(inverse-MEG-HEAT ${CMAKE_CURRENT_BINARY_DIR}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGHEAT.txt 1 1e-4 HEAT 10000 0.0)

ADD_TEST(inverse-MEG-MN ${CMAKE_CURRENT_BINARY_DIR}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGMN.txt 1 1e-4 MN 10000 0.0)

ADD_TEST(inverse-MEG-TV ${CMAKE_CURRENT_BINARY_DIR}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGTV.txt 1 1e-4 TV 10000 0.0)


############ TEST DIPOLE FORWARD RESULTS (Regression test) ##############

# assemble -rhsPOINT geometry.geom conductivity.cond dipoles.dip rhsmat.bin
ADD_TEST(assemble-rhsPOINT ${CMAKE_CURRENT_BINARY_DIR}/assemble -rhsPOINT ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.dip ${OpenMEEG_BINARY_DIR}/tests/rhspointmatrix.bin)

ADD_TEST(assemble-x2EEGrespPoint ${CMAKE_CURRENT_BINARY_DIR}/assemble -xToEEGresponse ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1_scalp.patches ${OpenMEEG_BINARY_DIR}/tests/x2eegpointmatrix.bin)

# assemble -sToMEGresponseContrib_point dipoles.dip squidscoord.squids sToMEGmat.bin
ADD_TEST(assemble-sToMEGrespPoint ${CMAKE_CURRENT_BINARY_DIR}/assemble -sToMEGresponseContrib_point ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.dip ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/s2megpointmatrix.bin)

ADD_TEST(gain-EEGPoint ${CMAKE_CURRENT_BINARY_DIR}/gain -EEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhspointmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2eegpointmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/eegpointgainmatrix.bin)

ADD_TEST(gain-MEGPoint ${CMAKE_CURRENT_BINARY_DIR}/gain -MEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhspointmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/s2megpointmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/megpointgainmatrix.bin)

# forward -EEG gainmatrix.bin dipoleActivation.src estimatedeegdata.txt noiselevel
ADD_TEST(forward-EEG-dipoles ${CMAKE_CURRENT_BINARY_DIR}/forward -EEG ${OpenMEEG_BINARY_DIR}/tests/eegpointgainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1_dip.src ${OpenMEEG_BINARY_DIR}/tests/estimatedEEGpoint.txt 0.0)

ADD_TEST(forward-MEG-dipoles ${CMAKE_CURRENT_BINARY_DIR}/forward -MEG ${OpenMEEG_BINARY_DIR}/tests/megpointgainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1_dip.src ${OpenMEEG_BINARY_DIR}/tests/estimatedMEGpoint.txt 0.0)

ENDIF(BUILD_TEST)
