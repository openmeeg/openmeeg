############ IO TEST ##############
ADD_EXECUTABLE(load_geo load_geo.cpp)
TARGET_LINK_LIBRARIES (load_geo ${OPENMEEG_OTHER_LIBRARIES} OpenMEEG)
ADD_TEST(load_geo ${CMAKE_CURRENT_BINARY_DIR}/load_geo ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.cond)

ADD_EXECUTABLE(test_sensors test_sensors.cpp)
TARGET_LINK_LIBRARIES (test_sensors OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(compare_matrix compare_matrix.cpp)
TARGET_LINK_LIBRARIES (compare_matrix ${OPENMEEG_OTHER_LIBRARIES})

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

SET(SUBJECT "Head1")

#------ LHS -----#
ADD_TEST(compareLHS-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                    ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.lhs
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.lhs
                    -sym)
#------ LHS INV -----#
ADD_TEST(compareLHSINV-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                       ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.lhs_inv
                       ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.lhs_inv
                       -sym)
#------ RHS -----#
ADD_TEST(compareRHS-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                    ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.rhs
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.rhs -full)

############ TEST EEG RESULTS ##############

#------ v2eeg -----#
ADD_TEST(comparev2eeg-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix 
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.v2eeg
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.v2eeg -if old_binary -sparse)
#------ EEG GAIN -----#
ADD_TEST(compareHeeg-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.Heeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.Heeg -full)
#------ EEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareEEGEST-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.est_eeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.est_eeg -full)
#------ EEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareEEG-HEAT-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                           ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-heat.est_src
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-heat.est_src -full -eps 0.02)
#------ EEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareEEG-MN-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-mn.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-mn.est_src -full)
#------ EEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareEEG-TV-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-tv.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-tv.est_src -full -eps 0.0001)

############ TEST MEG RESULTS ##############

#------ v2meg -----#
ADD_TEST(comparev2meg-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.v2meg
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.v2meg -full)
#------ S2MEG -----#
ADD_TEST(compareS2MEG-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.s2meg
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.s2meg -full)
#------ MEG GAIN -----#
ADD_TEST(compareMEGGAIN-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.Hmeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.Hmeg -full)
#------ MEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareMEGEST-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.est_meg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.est_meg -full)
#------ MEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareMEG-HEAT-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                           ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-heat.est_src
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-heat.est_src -full -eps 0.02)
#------ MEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareMEG-MN-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-mn.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-mn.est_src -full)
#------ MEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareMEG-TV-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-tv.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-tv.est_src -full)

############ TEST EEG RESULTS ON DIPOLES ##############

FOREACH ( DIP 1 2 3 4 5 )

    SET(EPSILON 0.15)

    # Compare EEG result with analytical solution obtained with Matlab
    ADD_TEST(compareEEGEST-point-Head1-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-point.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGEST-point-Head2-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-point.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head2_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    IF ( TEST_HEAD3 )
        ADD_TEST(compareEEGEST-point-Head3-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-point.est_eeg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head3_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP} -full)
    ENDIF ( TEST_HEAD3 )

ENDFOREACH ( DIP 1 2 3 4 5 )

FOREACH ( DIP 1 2 3 4 5)

    SET(EPSILON 0.15)

    # Compare MEG result with analytical solution obtained with Matlab
    ADD_TEST(compareMEGEST-point-Head1-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-point.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareMEGEST-point-Head2-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-point.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    IF ( TEST_HEAD3 )
        ADD_TEST(compareMEGEST-point-Head3-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-point.est_meg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP} -full)
    ENDIF ( TEST_HEAD3 )

ENDFOREACH ( DIP 1 2 3 4 5)

# Compare MEG result with analytical solution obtained with Matlab for dipole 6 that should be radial
SET(DIP 6)
SET(EPSILON 0.001)

ADD_TEST(compareMEGEST-point-Head1-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/Head1-point.est_meg
                         ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                         -eps ${EPSILON} -col ${DIP} -full)

ADD_TEST(compareMEGEST-point-Head2-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/Head2-point.est_meg
                         ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                         -eps ${EPSILON} -col ${DIP} -full)

IF ( TEST_HEAD3 )
    ADD_TEST(compareMEGEST-point-Head3-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head3-point.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -eps ${EPSILON} -col ${DIP} -full)
ENDIF ( TEST_HEAD3 )

# set tests that are expected to fail :
SET_TESTS_PROPERTIES(compareEEGEST-point-Head1-d4 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGEST-point-Head1-d5 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareMEGEST-point-Head1-d6 PROPERTIES WILL_FAIL TRUE)

