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
ADD_TEST(assemble-LHS ${EXECUTABLE_OUTPUT_PATH}/assemble -LHS ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin)

ADD_TEST(minverser-lhsmatrix ${EXECUTABLE_OUTPUT_PATH}/minverser ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin)

ADD_TEST(assemble-RHS ${EXECUTABLE_OUTPUT_PATH}/assemble -RHS ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin)

ADD_TEST(presmooth ${EXECUTABLE_OUTPUT_PATH}/presmooth ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin)


############ EEG TEST ##############
ADD_TEST(assemble-x2EEGresp ${EXECUTABLE_OUTPUT_PATH}/assemble -xToEEGresponse ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.patches ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin)

ADD_TEST(gain-EEG ${EXECUTABLE_OUTPUT_PATH}/gain -EEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin)

#------ FORWARD -----#
ADD_TEST(forward-EEG ${EXECUTABLE_OUTPUT_PATH}/forward -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.src ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt 0.0)

#------ INVERSE -----#
ADD_TEST(inverse-EEG-HEAT ${EXECUTABLE_OUTPUT_PATH}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGHEAT.txt 1 1e-4 HEAT 10000 0.0)

ADD_TEST(inverse-EEG-MN ${EXECUTABLE_OUTPUT_PATH}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGMN.txt 1 1e-4 MN 10000 0.0)

ADD_TEST(inverse-EEG-TV ${EXECUTABLE_OUTPUT_PATH}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGTV.txt 1 1e-4 TV 10000 0.0)


############ MEG TEST ##############
ADD_TEST(assemble-x2MEGresp ${EXECUTABLE_OUTPUT_PATH}/assemble -xToMEGresponseContrib ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin)

ADD_TEST(assemble-s2MEGresp ${EXECUTABLE_OUTPUT_PATH}/assemble -sToMEGresponseContrib ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin)

ADD_TEST(gain-MEG ${EXECUTABLE_OUTPUT_PATH}/gain -MEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin)

#------ FORWARD -----#
ADD_TEST(forward-MEG ${EXECUTABLE_OUTPUT_PATH}/forward -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.src ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt 0.0)

#------ INVERSE -----#
ADD_TEST(inverse-MEG-HEAT ${EXECUTABLE_OUTPUT_PATH}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGHEAT.txt 1 1e-4 HEAT 10000 0.0)

ADD_TEST(inverse-MEG-MN ${EXECUTABLE_OUTPUT_PATH}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGMN.txt 1 1e-4 MN 10000 0.0)

ADD_TEST(inverse-MEG-TV ${EXECUTABLE_OUTPUT_PATH}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGTV.txt 1 1e-4 TV 10000 0.0)

ENDIF(BUILD_TEST)
