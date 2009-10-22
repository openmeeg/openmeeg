INCLUDE(${OpenMEEG_SOURCE_DIR}/macros/ComparisonTest.cmake)

SET_FILE_PROPERTIES(Suffix
    "HeadMat .hm" "HeadMatInv .hm_inv" "SurfSourceMat .ssm" "v2eeg .h2em" "SurfGainEEG .sgem"
    "EEGEST .est_eeg" "EEG-HEAT -eeg-heat.est_src" "EEG-MN -eeg-mn.est_src" "EEG-TV -eeg-tv.est_src"
    "Head2MEGMat .h2mm" "SurfSource2MEGMat .ss2mm" "SurfGainMEG .sgmm" "MEGEST .est_meg"
    "MEG-HEAT -meg-heat.est_src" "MEG-MN -meg-mn.est_src" "MEG-TV -meg-tv.est_src"
)

SET_FILE_PROPERTIES(CompareOptions
    "HeadMat -sym" "HeadMatInv -sym" "SurfSourceMat -full:-if1:binary:-if2:binary" "v2eeg -if:binary:-sparse"
    "SurfGainEEG -full" "EEGEST -full" "EEG-HEAT -full:-eps:0.02" "EEG-MN -full" "EEG-TV -full:-eps:0.0001"
    "Head2MEGMat -full" "SurfSource2MEGMat -full" "SurfGainMEG -full" "MEGEST -full"
    "MEG-HEAT -full:-eps:0.02" "MEG-MN -full" "MEG-TV -full"
)

############ IO TEST ##############
OPENMEEG_UNIT_TEST(load_geo
    SOURCES load_geo.cpp
    LIBRARIES OpenMEEG OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}"
    PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.cond)

NEW_EXECUTABLE(test_sensors test_sensors.cpp
               LIBRARIES OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")
NEW_EXECUTABLE(compare_matrix compare_matrix.cpp
               LIBRARIES OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}" ${LAPACK_LIBRARIES})
NEW_EXECUTABLE(om_validationEIT validationEIT.cpp ${OPEMEEG_HEADERS}
               LIBRARIES OpenMEEG OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}" ${LAPACK_LIBRARIES})

IF ( WIN32 )
    SET(COMPARE_MATRIX "${EXECUTABLE_OUTPUT_PATH}/compare_matrix")
    SET(VALIDATION_EIT "${EXECUTABLE_OUTPUT_PATH}/om_validationEIT")
ELSE ( WIN32 )
    SET(COMPARE_MATRIX "${CMAKE_CURRENT_BINARY_DIR}/compare_matrix")
    SET(VALIDATION_EIT "${CMAKE_CURRENT_BINARY_DIR}/om_validationEIT")
ENDIF ( WIN32 )

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

SET(COMPARISONS HeadMat HeadMatInv SurfSourceMat v2eeg SurfGainEEG EEGEST EEG-HEAT EEG-MN EEG-TV
                Head2MEGMat SurfSource2MEGMat SurfGainMEG MEGEST MEG-HEAT MEG-MN MEG-TV
)

FOREACH (COMPARISON ${COMPARISONS})
    OPENMEEG_COMPARISON_TEST(Head1 ${COMPARISON} compare_matrix)
ENDFOREACH()

############ TEST EEG RESULTS ON DIPOLES ##############

FOREACH ( DIP 1 2 3 4 5 )

    SET(EPSILON 0.15)

    # Compare EEG result with analytical solution obtained with Matlab
    ADD_TEST(compareEEGEST-dip-Head1-d${DIP}-rdm ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGEST-dip-Head2-d${DIP}-rdm ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head2_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGEST-dip-Head1-d${DIP}-mag ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head1_analytic.txt
                             -mag -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGEST-dip-Head2-d${DIP}-mag ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_eeg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head2_analytic.txt
                             -mag -eps ${EPSILON} -col ${DIP} -full)
    IF ( TEST_HEAD3 )
        ADD_TEST(compareEEGEST-dip-Head3-d${DIP}-rdm ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_eeg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head3_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP} -full)

        ADD_TEST(compareEEGEST-dip-Head3-d${DIP}-mag ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_eeg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_head3_analytic.txt
                                 -mag -eps ${EPSILON} -col ${DIP} -full)
    ENDIF ( TEST_HEAD3 )

ENDFOREACH ( DIP 1 2 3 4 5 )

FOREACH ( SUB 1 2 )

    # validationEIT geometry.geom conductivity.cond dipoles.dip source.dsm headmatinv.bin out.eit_qgradVj out.eit_diffVf
    ADD_TEST(EITvalidation-dipoles-Head${SUB}
        ${VALIDATION_EIT}
        ${OpenMEEG_SOURCE_DIR}/data/Models/Head${SUB}/Head${SUB}.geom
        ${OpenMEEG_SOURCE_DIR}/data/Models/Head${SUB}/Head${SUB}.cond
        ${OpenMEEG_SOURCE_DIR}/data/Computations/Head${SUB}/Head${SUB}.dip
        ${OpenMEEG_BINARY_DIR}/tests/Head${SUB}.dsm
        ${OpenMEEG_BINARY_DIR}/tests/Head${SUB}.hm_inv
        ${OpenMEEG_BINARY_DIR}/tests/Head${SUB}.eit_qgradVj
        ${OpenMEEG_BINARY_DIR}/tests/Head${SUB}.eit_diffVf)

ENDFOREACH ( SUB 1 2 )

FOREACH ( DIP 1 2 3 4 5 6 )

    SET(EPSILON 0.15)
    # Compare the potential results in a interior sphere of the Surf2Vol operator with analytical solution
    # obtained with Sphere (V.Hedou Modified)
    ADD_TEST(compareEEGinternal-dip-Head1-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_eeginternal
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_internal_analytic.txt
                             -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEEGinternal-dip-Head2-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_eeginternal
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_internal_analytic.txt
                             -eps ${EPSILON} -col ${DIP} -full)

    # Compare the q.gradVj to diffVf in order to validate the EIT problem
    ADD_TEST(compareEITvalidation-Head1-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1.eit_qgradVj
                             ${OpenMEEG_BINARY_DIR}/tests/Head1.eit_diffVf
                             -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareEITvalidation-Head2-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2.eit_qgradVj
                             ${OpenMEEG_BINARY_DIR}/tests/Head2.eit_diffVf
                             -eps ${EPSILON} -col ${DIP} -full)

ENDFOREACH ( DIP 1 2 3 4 5 6 )

FOREACH ( DIP 1 2 3 4 5)

    SET(EPSILON 0.15)

    # Compare MEG result with analytical solution obtained with Matlab
    ADD_TEST(compareMEGEST-dip-Head1-d${DIP}-rdm ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareMEGEST-dip-Head2-d${DIP}-rdm ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -rdm -eps ${EPSILON} -col ${DIP} -full)

    IF ( TEST_HEAD3 )
        ADD_TEST(compareMEGEST-dip-Head3-d${DIP}-rdm ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_meg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                                 -rdm -eps ${EPSILON} -col ${DIP} -full)
    ENDIF ( TEST_HEAD3 )

    # Compare MEG result with analytical solution obtained with Matlab
    ADD_TEST(compareMEGEST-dip-Head1-d${DIP}-mag ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head1-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -mag -eps ${EPSILON} -col ${DIP} -full)

    ADD_TEST(compareMEGEST-dip-Head2-d${DIP}-mag ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head2-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -mag -eps ${EPSILON} -col ${DIP} -full)

    IF ( TEST_HEAD3 )
        ADD_TEST(compareMEGEST-dip-Head3-d${DIP}-mag ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head3-dip.est_meg
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                                 -mag -eps ${EPSILON} -col ${DIP} -full)
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
SET_TESTS_PROPERTIES(compareEEGEST-dip-Head1-d1-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGEST-dip-Head1-d4-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGEST-dip-Head1-d5-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d2 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d3 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d4 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d5 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d6 PROPERTIES WILL_FAIL TRUE)
