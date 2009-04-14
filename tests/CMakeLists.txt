############ IO TEST ##############
ADD_EXECUTABLE(load_geo load_geo.cpp)
TARGET_LINK_LIBRARIES (load_geo ${OPENMEEG_OTHER_LIBRARIES} OpenMEEG)

ADD_EXECUTABLE(test_sensors test_sensors.cpp)
TARGET_LINK_LIBRARIES (test_sensors OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(compare_matrix compare_matrix.cpp)
TARGET_LINK_LIBRARIES (compare_matrix ${OPENMEEG_OTHER_LIBRARIES})

IF ( WIN32 )
    SET(LOAD_GEO "${EXECUTABLE_OUTPUT_PATH}/load_geo")
    SET(COMPARE_MATRIX "${EXECUTABLE_OUTPUT_PATH}/compare_matrix")
ELSE ( WIN32 )
    SET(LOAD_GEO "${CMAKE_CURRENT_BINARY_DIR}/load_geo")
    SET(COMPARE_MATRIX "${CMAKE_CURRENT_BINARY_DIR}/compare_matrix")
ENDIF ( WIN32 )

ADD_TEST(load_geo ${LOAD_GEO} ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.cond)

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

SET(SUBJECT "Head1")

#------ HeadMat -----#
ADD_TEST(compareHeadMat-${SUBJECT} ${COMPARE_MATRIX}
                    ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.hm
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.hm
                    -sym)
#------ HeadMat INV -----#
ADD_TEST(compareHeadMatINV-${SUBJECT} ${COMPARE_MATRIX}
                       ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.hm_inv
                       ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.hm_inv
                       -sym)
#------ SurfSourceMat -----#
ADD_TEST(compareSurfSourceMat-${SUBJECT} ${COMPARE_MATRIX}
                    ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.ssm
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.ssm -full -if1 binary -if2 binary)

############ TEST EEG RESULTS ##############

#------ v2eeg -----#
ADD_TEST(comparev2eeg-${SUBJECT} ${COMPARE_MATRIX} 
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.h2em
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.h2em -if binary -sparse)
#------ EEG GAIN -----#
ADD_TEST(compareSurfGainEEG-${SUBJECT} ${COMPARE_MATRIX}
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.sgem
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.sgem -full)
#------ EEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareEEGEST-${SUBJECT} ${COMPARE_MATRIX}
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.est_eeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.est_eeg -full)
#------ EEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareEEG-HEAT-${SUBJECT} ${COMPARE_MATRIX}
                           ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-heat.est_src
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-heat.est_src -full -eps 0.02)
#------ EEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareEEG-MN-${SUBJECT} ${COMPARE_MATRIX}
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-mn.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-mn.est_src -full)
#------ EEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareEEG-TV-${SUBJECT} ${COMPARE_MATRIX}
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-tv.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-tv.est_src -full -eps 0.0001)

############ TEST MEG RESULTS ##############

#------ Head2MEGMat -----#
ADD_TEST(compareHead2MEGMat-${SUBJECT} ${COMPARE_MATRIX}
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.h2mm
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.h2mm -full)
#------ SurfSource2MEGMat -----#
ADD_TEST(compareSurfSource2MEGMat-${SUBJECT} ${COMPARE_MATRIX}
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.ss2mm
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.ss2mm -full)
#------ MEG GAIN -----#
ADD_TEST(compareSurfGainMEG-${SUBJECT} ${COMPARE_MATRIX}
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.sgmm
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.sgmm -full)
#------ MEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareMEGEST-${SUBJECT} ${COMPARE_MATRIX}
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.est_meg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.est_meg -full)
#------ MEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareMEG-HEAT-${SUBJECT} ${COMPARE_MATRIX}
                           ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-heat.est_src
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-heat.est_src -full -eps 0.02)
#------ MEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareMEG-MN-${SUBJECT} ${COMPARE_MATRIX}
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-mn.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-mn.est_src -full)
#------ MEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareMEG-TV-${SUBJECT} ${COMPARE_MATRIX}
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-tv.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-tv.est_src -full)

############ TEST EEG RESULTS ON DIPOLES ##############

FOREACH ( DIP 1 2 3 4 5 )

    SET(EPSILON 0.15)

    # Compare EEG result with analytical solution obtained with Matlab
    ADD_TEST(compareEEGEST-dip-Head1-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGEST-dip-Head2-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head2_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    IF ( TEST_HEAD3 )
        ADD_TEST(compareEEGEST-dip-Head3-d${DIP} ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_eeg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head3_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP} -full)
    ENDIF ( TEST_HEAD3 )

ENDFOREACH ( DIP 1 2 3 4 5 )

FOREACH ( DIP 1 2 3 4 5 6 )

    SET(EPSILON 0.15)
    # Compare the potential results in a interior sphere of the Surf2Vol operator with analytical solution
    # obtained with Sphere (V.Hedou Modified)
    ADD_TEST(compareEEGINTERNAL-dip-Head1-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_eeginternal
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_internal_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGINTERNAL-dip-Head2-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_eeginternal
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_internal_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)
ENDFOREACH ( DIP 1 2 3 4 5 6 )


FOREACH ( DIP 1 2 3 4 5)

    SET(EPSILON 0.15)

    # Compare MEG result with analytical solution obtained with Matlab
    ADD_TEST(compareMEGEST-dip-Head1-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareMEGEST-dip-Head2-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    IF ( TEST_HEAD3 )
        ADD_TEST(compareMEGEST-dip-Head3-d${DIP} ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_meg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP} -full)
    ENDIF ( TEST_HEAD3 )

ENDFOREACH ( DIP 1 2 3 4 5)

# Compare MEG result with analytical solution obtained with Matlab for dipole 6 that should be radial
SET(DIP 6)
SET(EPSILON 0.001)

ADD_TEST(compareMEGEST-dip-Head1-d${DIP} ${COMPARE_MATRIX}
                         ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_meg
                         ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                         -eps ${EPSILON} -col ${DIP} -full)

ADD_TEST(compareMEGEST-dip-Head2-d${DIP} ${COMPARE_MATRIX}
                         ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_meg
                         ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                         -eps ${EPSILON} -col ${DIP} -full)

IF ( TEST_HEAD3 )
    ADD_TEST(compareMEGEST-dip-Head3-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -eps ${EPSILON} -col ${DIP} -full)
ENDIF ( TEST_HEAD3 )

# set tests that are expected to fail :
SET_TESTS_PROPERTIES(compareEEGEST-dip-Head1-d4 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGEST-dip-Head1-d5 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareMEGEST-dip-Head1-d6 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGINTERNAL-dip-Head1-d1 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGINTERNAL-dip-Head1-d2 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGINTERNAL-dip-Head1-d3 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGINTERNAL-dip-Head1-d4 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGINTERNAL-dip-Head1-d5 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGINTERNAL-dip-Head1-d6 PROPERTIES WILL_FAIL TRUE)

