############ IO TEST ##############
ADD_EXECUTABLE(load_geo load_geo.cpp)
TARGET_LINK_LIBRARIES (load_geo ${OPENMEEG_OTHER_LIBRARIES} OpenMEEG)
ADD_TEST(load_mesh_tri ${EXECUTABLE_OUTPUT_PATH}/load_geo ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri)

ADD_EXECUTABLE(sensors sensors.cpp)
TARGET_LINK_LIBRARIES (sensors ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(compare_matrix compare_matrix.cpp)
TARGET_LINK_LIBRARIES (compare_matrix ${OPENMEEG_OTHER_LIBRARIES})

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

SET(SUBJECT "Head1")

#------ LHS -----#
ADD_TEST(compareLHS-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                    ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.lhs
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.lhs
                    -sym -bin)
#------ LHS INV -----#
ADD_TEST(compareLHSINV-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                       ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.lhs_inv
                       ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.lhs_inv
                       -sym -bin)
#------ RHS -----#
ADD_TEST(compareRHS-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                    ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.rhs
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.rhs -bin)

############ TEST EEG RESULTS ##############

#------ X2EEG -----#
ADD_TEST(compareX2EEG-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.x2eeg
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.x2eeg -bin)
#------ EEG GAIN -----#
ADD_TEST(compareHeeg-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.Heeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.Heeg -bin)
#------ EEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareEEGEST-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.est_eeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.est_eeg -txt)
#------ EEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareEEG-HEAT-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                           ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-heat.est_src
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-heat.est_src -txt)
#------ EEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareEEG-MN-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-mn.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-mn.est_src -txt)
#------ EEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareEEG-TV-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-eeg-tv.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-eeg-tv.est_src -txt)

############ TEST MEG RESULTS ##############

#------ X2MEG -----#
ADD_TEST(compareX2MEG ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.x2meg
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.x2meg -bin)
#------ S2MEG -----#
ADD_TEST(compareS2MEG ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.s2meg
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.s2meg -bin)
#------ MEG GAIN -----#
ADD_TEST(compareMEGGAIN ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.Hmeg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.Hmeg -bin)
#------ MEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareMEGEST ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}.est_meg
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}.est_meg -txt)
#------ MEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareMEG-HEAT-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                           ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-heat.est_src
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-heat.est_src -txt)
#------ MEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareMEG-MN-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-mn.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-mn.est_src -txt)
#------ MEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareMEG-TV-${SUBJECT} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/${SUBJECT}-meg-tv.est_src
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/${SUBJECT}-meg-tv.est_src -txt)

############ TEST EEG RESULTS ON DIPOLES ##############

FOREACH ( DIP 1 2 3 4 5 )

    SET(EPSILON 0.15)

    # Compare EEG result with analytical solution obtained with Matlab
    ADD_TEST(compareEEGEST-point-Head1-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-point.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP})

    ADD_TEST(compareEEGEST-point-Head2-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-point.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head2_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP})

    IF ( TEST_HEAD3 )
        ADD_TEST(compareEEGEST-point-Head3-d${DIP} ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-point.est_eeg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head3_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP})
    ENDIF ( TEST_HEAD3 )

ENDFOREACH ( DIP 1 2 3 4 5 )
