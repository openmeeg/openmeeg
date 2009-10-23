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
    SET(BASE_FILE_NAME Head1${Suffix_${COMPARISON}})
    OPENMEEG_COMPARISON_TEST(Head1 ${COMPARISON} "" compare_matrix ${BASE_FILE_NAME} initialTest/${BASE_FILE_NAME})
ENDFOREACH()

############ TEST EEG RESULTS ON DIPOLES ##############

IF (TEST_HEAD3)
    SET(HEAD3 3)
ENDIF()

SET(EPSILON 0.15)
FOREACH(DIP 1 2 3 4 5)
    FOREACH(HEADNUM 1 2 ${HEAD3})
        FOREACH(COMP mag rdm)
            SET(BASE_FILE_NAME "Head${HEADNUM}-dip.est_eeg")
            # Compare EEG result with analytical solution obtained with Matlab
            OPENMEEG_COMPARISON_TEST(Head${HEADNUM}-dip EEGEST-dip "${DIP}-${COMP}"
                                     compare_matrix ${BASE_FILE_NAME} analytic/eeg_head${HEADNUM}_analytic.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full)
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

FOREACH(HEADNUM 1 2)
    # validationEIT geometry.geom conductivity.cond dipoles.dip source.dsm headmatinv.bin out.eit_qgradVj out.eit_diffVf
    ADD_TEST(EITvalidation-dipoles-Head${HEADNUM}
        ${VALIDATION_EIT}
        ${OpenMEEG_SOURCE_DIR}/data/Models/Head${HEADNUM}/Head${HEADNUM}.geom
        ${OpenMEEG_SOURCE_DIR}/data/Models/Head${HEADNUM}/Head${HEADNUM}.cond
        ${OpenMEEG_SOURCE_DIR}/data/Computations/Head${HEADNUM}/Head${HEADNUM}.dip
        ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}.dsm
        ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}.hm_inv
        ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}.eit_qgradVj
        ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}.eit_diffVf)
ENDFOREACH()

FOREACH(DIP 1 2 3 4 5 6)
    FOREACH(HEADNUM 1 2 ${HEAD3})
        # Compare the potential results in a interior sphere of the Surf2Vol operator with analytical solution
        # obtained with Sphere (V.Hedou Modified)
        ADD_TEST(compareEEGinternal-dip-Head${HEADNUM}-d${DIP} ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}-dip.est_eeginternal
                                 ${OpenMEEG_SOURCE_DIR}/tests/analytic/eeg_internal_analytic.txt
                                 -eps ${EPSILON} -col ${DIP} -full)

        # Compare the q.gradVj to diffVf in order to validate the EIT problem
        ADD_TEST(compareEITvalidation-Head${HEADNUM}-d${DIP} ${COMPARE_MATRIX}
                                 ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}.eit_qgradVj
                                 ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}.eit_diffVf
                                 -eps ${EPSILON} -col ${DIP} -full)
    ENDFOREACH()
ENDFOREACH()

FOREACH(DIP 1 2 3 4 5)
    FOREACH(HEADNUM 1 2 ${HEAD3})
        FOREACH(COMP mag rdm)
            # Compare MEG result with analytical solution obtained with Matlab
            ADD_TEST(compareMEGEST-dip-Head${HEADNUM}-d${DIP}-${COMP} ${COMPARE_MATRIX}
                                     ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}-dip.est_meg
                                     ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                                     -${COMP} -eps ${EPSILON} -col ${DIP} -full)
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

# Compare MEG result with analytical solution obtained with Matlab for dipole 6 that should be radial
SET(DIP 6)
SET(EPSILON 0.001)

FOREACH(HEADNUM 1 2 ${HEAD3})
    ADD_TEST(compareMEGEST-dip-Head${HEADNUM}-d${DIP} ${COMPARE_MATRIX}
                             ${OpenMEEG_BINARY_DIR}/tests/Head${HEADNUM}-dip.est_meg
                             ${OpenMEEG_SOURCE_DIR}/tests/analytic/meg_head1_analytic.txt
                             -eps ${EPSILON} -col ${DIP} -full)
ENDFOREACH()

# set tests that are expected to fail :
SET_TESTS_PROPERTIES(compare-EEGEST-dip-Head1-dip1-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compare-EEGEST-dip-Head1-dip4-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compare-EEGEST-dip-Head1-dip5-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d2 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d3 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d4 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d5 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(compareEEGinternal-dip-Head1-d6 PROPERTIES WILL_FAIL TRUE)
