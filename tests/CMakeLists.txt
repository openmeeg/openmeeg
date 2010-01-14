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

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

SET(COMPARISONS HeadMat HeadMatInv SurfSourceMat v2eeg SurfGainEEG EEGEST EEG-HEAT EEG-MN EEG-TV
                Head2MEGMat SurfSource2MEGMat SurfGainMEG MEGEST MEG-HEAT MEG-MN MEG-TV
)

FOREACH (COMPARISON ${COMPARISONS})
    SET(BASE_FILE_NAME Head1${Suffix_${COMPARISON}})
    SET(TEST_NAME ${COMPARISON}-Head1)
    OPENMEEG_COMPARISON_TEST(${COMPARISON}-Head1 ${BASE_FILE_NAME} initialTest/${BASE_FILE_NAME} ${CompareOptions_${COMPARISON}})
ENDFOREACH()

############ TEST EEG RESULTS ON DIPOLES ##############

IF (TEST_HEAD3)
    SET(HEAD3 3)
ENDIF()

SET(EPSILON 0.15)
FOREACH(DIP 1 2 3 4 5)
    FOREACH(HEADNUM 1 2 ${HEAD3})
        FOREACH(COMP mag rdm)
            SET(HEAD "Head${HEADNUM}")
            SET(BASE_FILE_NAME "${HEAD}-dip.est_eeg")
            # Compare EEG result with analytical solution obtained with Matlab
            OPENMEEG_COMPARISON_TEST("EEGEST-dip-${HEAD}-dip${DIP}-${COMP}"
                                     ${BASE_FILE_NAME} analytic/eeg_head${HEADNUM}_analytic.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full)
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()
# set tests that are expected to fail :
SET_TESTS_PROPERTIES(cmp-EEGEST-dip-Head1-dip1-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-EEGEST-dip-Head1-dip4-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-EEGEST-dip-Head1-dip5-rdm PROPERTIES WILL_FAIL TRUE)

IF (WIN32)
    SET(VALIDATION_EIT "${EXECUTABLE_OUTPUT_PATH}/om_validationEIT")
ELSE()
    SET(VALIDATION_EIT "${CMAKE_CURRENT_BINARY_DIR}/om_validationEIT")
ENDIF()

FOREACH(HEADNUM 1 2)
    SET(HEAD "Head${HEADNUM}")
    SET(DATADIR ${OpenMEEG_SOURCE_DIR}/data)
    SET(MODELPREFIX ${DATADIR}/Models/${HEAD}/${HEAD})
    SET(TESTPREFIX ${OpenMEEG_BINARY_DIR}/tests/${HEAD})
    # validationEIT geometry.geom conductivity.cond dipoles.dip source.dsm headmatinv.bin out.eit_qgradVj out.eit_diffVf
    ADD_TEST(EITvalidation-dipoles-${HEAD}
        ${VALIDATION_EIT} ${MODELPREFIX}.geom ${MODELPREFIX}.cond
        ${DATADIR}/Computations/${HEAD}/${HEAD}.dip
        ${TESTPREFIX}.dsm ${TESTPREFIX}.hm_inv ${TESTPREFIX}.eit_qgradVj ${TESTPREFIX}.eit_diffVf)
ENDFOREACH()

FOREACH(DIP 1 2 3 4 5 6)
    FOREACH(HEADNUM 1 2)
        SET(HEAD "Head${HEADNUM}")
        # Compare the potential results in a interior sphere of the Surf2Vol operator with analytical solution
        # obtained with Sphere (V.Hedou Modified)
        OPENMEEG_COMPARISON_TEST("EEGinternal-dip-${HEAD}-dip${DIP}"
                                 ${HEAD}-dip.est_eeginternal analytic/eeg_internal_analytic.txt -eps ${EPSILON} -col ${DIP} -full)

        # Compare the q.gradVj to diffVf in order to validate the EIT problem
        OPENMEEG_COMPARISON_TEST("EIT-${HEAD}-dip${DIP}"
                                 ${HEAD}.eit_qgradVj ${OpenMEEG_BINARY_DIR}/tests/${HEAD}.eit_diffVf -eps ${EPSILON} -col ${DIP} -full)
    ENDFOREACH()
ENDFOREACH()
# set tests that are expected to fail :
SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head1-dip2 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head1-dip3 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head1-dip4 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head1-dip5 PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head1-dip6 PROPERTIES WILL_FAIL TRUE)

FOREACH(SENSORORIENT "" "-tangential" "-noradial")
    FOREACH(HEADNUM 1 2 ${HEAD3})
        SET(HEAD "Head${HEADNUM}")
        FOREACH(DIP 1 2 3 4 5)
            FOREACH(COMP mag rdm)
                # Compare MEG result with analytical solution obtained with Matlab
                OPENMEEG_COMPARISON_TEST("MEGEST-dip-${HEAD}-dip${DIP}${SENSORORIENT}-${COMP}"
                    ${HEAD}-dip${SENSORORIENT}.est_meg analytic/meg_analytic${SENSORORIENT}.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full)
            ENDFOREACH()
        ENDFOREACH()
        OPENMEEG_COMPARISON_TEST("MEGEST-dip${SENSORORIENT}-${HEAD}-dip6"
            ${HEAD}-dip${SENSORORIENT}.est_meg analytic/meg_analytic${SENSORORIENT}.txt -eps 0.001 -col 6 -full)
    ENDFOREACH()
ENDFOREACH()
# set tests that are expected to fail :
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip2-tangential-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip2-tangential-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip3-tangential-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip3-tangential-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip4-tangential-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip4-tangential-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip5-tangential-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head2-dip4-tangential-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head2-dip4-tangential-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head2-dip5-tangential-mag PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head2-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip4-noradial-mag PROPERTIES   WILL_FAIL TRUE)
SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head1-dip5-noradial-mag PROPERTIES   WILL_FAIL TRUE)
