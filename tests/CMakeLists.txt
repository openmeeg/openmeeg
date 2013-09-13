INCLUDE(ComparisonTest)

SET_FILE_PROPERTIES(Suffix
    "HM .hm" "HMInv .hm_inv" "DSM .dsm" "SSM .ssm" "v2eeg .h2em" "SurfGainEEG .sgem"
    "ESTEEG .est_eeg" "EEGadjointEST .est_eegadjoint" "MEGadjointEST .est_megadjoint"
    "Head2MEGMat .h2mm" "SS2MM .ss2mm" "SurfGainMEG .sgmm" "ESTMEG .est_meg"
)

SET_FILE_PROPERTIES(CompareOptions
    "HM -sym" "HMInv -sym" "DSM -full:-if1:binary:-if2:binary" "SSM -full:-if1:binary:-if2:binary" "v2eeg -if:binary:-sparse"
    "SurfGainEEG -full" "ESTEEG -full" "Head2MEGMat -full" "SS2MM -full" "SurfGainMEG -full" "ESTMEG -full"
)

############ IO TEST ##############
OPENMEEG_UNIT_TEST(load_geo
    SOURCES load_geo.cpp
    LIBRARIES OpenMEEG OpenMEEGMaths
    PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.cond)

OPENMEEG_UNIT_TEST(test_mesh_ios
    SOURCES test_mesh_ios.cpp
    LIBRARIES OpenMEEG OpenMEEGMaths
    PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Models/Head1/Head1.tri)

NEW_EXECUTABLE(test_sensors test_sensors.cpp
               LIBRARIES OpenMEEG)
NEW_EXECUTABLE(compare_matrix compare_matrix.cpp
               LIBRARIES OpenMEEGMaths ${LAPACK_LIBRARIES})
NEW_EXECUTABLE(om_validationEIT validationEIT.cpp ${OPEMEEG_HEADERS}
               LIBRARIES OpenMEEG OpenMEEGMaths ${LAPACK_LIBRARIES})

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

SET(COMPARISONS HM HMInv SSM DSM v2eeg SurfGainEEG ESTEEG
                Head2MEGMat SS2MM SurfGainMEG ESTMEG 
)

FOREACH (COMPARISON ${COMPARISONS})
    SET(BASE_FILE_NAME Head1${Suffix_${COMPARISON}})
    SET(TEST_NAME ${COMPARISON}-Head1)
    OPENMEEG_COMPARISON_TEST(${COMPARISON}-Head1
        ${BASE_FILE_NAME} initialTest/${BASE_FILE_NAME} ${CompareOptions_${COMPARISON}})
ENDFOREACH()

############ TEST EEG RESULTS ON DIPOLES ##############

# defining variables for those who do not use VTK
IF (USE_VTK) 
    SET(NN "NNa" "NNb" "NNc")
    SET(NNa "NNa")
    SET(NNb "NNb")
    SET(NNc "NNc")
    SET(NNa1 "NNa1")
    SET(NNb1 "NNb1")
    SET(NNc1 "NNc1")
    SET(NNa2 "NNa2")
    SET(NNb2 "NNb2")
    SET(NNc2 "NNc2")
ENDIF()

IF (TEST_HEAD3)
    SET(HEAD3 3)
ENDIF()

SET(EPSILON 0.15)
FOREACH(DIP 1 2 3 4 5)
    FOREACH(HEADGEO "" ${NN})
        FOREACH(HEADNUM 1 2 ${HEAD3})
            FOREACH(COMP mag rdm)
                SET(HEAD "Head${HEADGEO}${HEADNUM}")
                FOREACH(ADJOINT "" adjoint adjoint2)
                    SET(BASE_FILE_NAME "${HEAD}-dip.est_eeg${ADJOINT}")
                    # Compare EEG result with analytical solution obtained with Matlab
                    OPENMEEG_COMPARISON_TEST("EEG${ADJOINT}EST-dip-${HEAD}-dip${DIP}-${COMP}"
                        ${BASE_FILE_NAME} analytic/eeg_head${HEADNUM}_analytic.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full
                        DEPENDS EEG${ADJOINT}-dipoles-${HEAD})
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()
IF (TEST_HEAD3)
    FOREACH(DIP 1 2)
        FOREACH(COMP mag rdm)
            SET(BASE_FILE_NAME "Head3-dip-skullscalp.est_eeg")
            # Compare EEG result with a numerical solution obtained with the tetrahedral FEM
            OPENMEEG_COMPARISON_TEST("EEGEST-dipSkullScalp-Head3-dip${DIP}-${COMP}"
                ${BASE_FILE_NAME} analytic/eeg_head3-skullscalp_FEM.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full
                DEPENDS EEG-dipolesSkullScalp-Head3)
        ENDFOREACH()
    ENDFOREACH()
ENDIF()

# set tests that are expected to fail :
FOREACH(ADJOINT "" "adjoint" "adjoint2")
    FOREACH(HEADGEO ${NNc1})
        FOREACH(DIP 1 2 3 4 5)
            SET_TESTS_PROPERTIES(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip${DIP}-mag PROPERTIES WILL_FAIL TRUE) # all cmp-EEG-mag NNc1 tests fail...
        ENDFOREACH()
    ENDFOREACH()
    FOREACH(HEADGEO "1" ${NNa1} ${NNb1})
        SET_TESTS_PROPERTIES(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip4-rdm PROPERTIES WILL_FAIL TRUE)
        SET_TESTS_PROPERTIES(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip5-rdm PROPERTIES WILL_FAIL TRUE)
    ENDFOREACH()
ENDFOREACH()

IF (WIN32)
    SET(VALIDATION_EIT "${EXECUTABLE_OUTPUT_PATH}/om_validationEIT")
ELSE()
    SET(VALIDATION_EIT "${CMAKE_CURRENT_BINARY_DIR}/om_validationEIT")
ENDIF()

FOREACH(HEADGEO "" ${NN})
    FOREACH(HEADNUM 1 2 ${HEAD3})
        SET(HEAD "Head${HEADGEO}${HEADNUM}")
        SET(DATADIR ${OpenMEEG_SOURCE_DIR}/data)
        SET(MODELPREFIX ${DATADIR}/Models/${HEAD}/${HEAD})
        SET(TESTPREFIX ${OpenMEEG_BINARY_DIR}/tests/${HEAD})
        # validationEIT geometry.geom conductivity.cond dipoles.dip source.dsm headmatinv.bin out.eit_qgradVj out.eit_diffVf
        OPENMEEG_TEST(
            EITvalidation-dipoles-${HEAD}
            ${VALIDATION_EIT} ${MODELPREFIX}.geom ${MODELPREFIX}.cond
            ${DATADIR}/Computations/${HEAD}/${HEAD}.dip
            ${TESTPREFIX}.dsm ${TESTPREFIX}.hm_inv ${TESTPREFIX}.eit_qgradVj ${TESTPREFIX}.eit_diffVf
            DEPENDS HMINV-${HEAD} DSM-${HEAD})
    ENDFOREACH()
ENDFOREACH()

FOREACH(DIP 1 2 3 4 5)
    FOREACH(HEADGEO "" ${NN})
        FOREACH(HEADNUM 1 2 ${HEAD3})
            SET(HEAD "Head${HEADGEO}${HEADNUM}")
            # Compare the potential results in a interior sphere of the Surf2Vol operator with analytical solution
            # obtained with Sphere (V.Hedou Modified)
            OPENMEEG_COMPARISON_TEST("EEGinternal-dip-${HEAD}-dip${DIP}"
                ${HEAD}-dip-internal.est_eeg analytic/eeg_internal_analytic.txt -eps ${EPSILON} -col ${DIP} -full
                DEPENDS InternalPot-dipoles-${HEAD})

            # Compare the q.gradVj to diffVf in order to validate the EIT problem
            OPENMEEG_COMPARISON_TEST("EIT-${HEAD}-dip${DIP}"
                ${HEAD}.eit_qgradVj ${OpenMEEG_BINARY_DIR}/tests/${HEAD}.eit_diffVf -eps ${EPSILON} -col ${DIP} -full
                DEPENDS EITvalidation-dipoles-${HEAD})
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

# set tests that are expected to fail :
FOREACH(HEADGEO "1" ${NNa1} ${NNb1})
    SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head${HEADGEO}-dip5 PROPERTIES WILL_FAIL TRUE)
ENDFOREACH()
FOREACH(HEADGEO ${NNa1})
    SET_TESTS_PROPERTIES(cmp-EIT-Head${HEADGEO}-dip1 PROPERTIES WILL_FAIL TRUE) # all cmp-EIT-HeadNNa1 tests fail...
    SET_TESTS_PROPERTIES(cmp-EIT-Head${HEADGEO}-dip2 PROPERTIES WILL_FAIL TRUE)
ENDFOREACH()
FOREACH(HEADGEO ${NNa1} ${NNa2})
    FOREACH(DIP 3 4 5)
        SET_TESTS_PROPERTIES(cmp-EIT-Head${HEADGEO}-dip${DIP} PROPERTIES WILL_FAIL TRUE)
    ENDFOREACH()
ENDFOREACH()

FOREACH(SENSORORIENT "" "-tangential" "-noradial")
    FOREACH(ADJOINT "" adjoint adjoint2)
        IF (NOT(${ADJOINT} STREQUAL "adjoint" OR ${ADJOINT} STREQUAL "adjoint2") OR (SENSORORIENT STREQUAL ""))
            FOREACH(HEADGEO "" ${NN})
                FOREACH(HEADNUM 1 2 ${HEAD3})
                    SET(HEAD "Head${HEADGEO}${HEADNUM}")
                    FOREACH(DIP 1 2 3 4 5)
                        FOREACH(COMP mag rdm)
                            # Compare MEG result with analytical solution obtained with Matlab
                            OPENMEEG_COMPARISON_TEST("MEG${ADJOINT}EST-dip-${HEAD}-dip${DIP}${SENSORORIENT}-${COMP}"
                                ${HEAD}-dip${SENSORORIENT}.est_meg${ADJOINT} analytic/meg_analytic${SENSORORIENT}.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full
                                DEPENDS MEG${ADJOINT}-dipoles-${HEAD}${SENSORORIENT})
                        ENDFOREACH()
                    ENDFOREACH()
                    OPENMEEG_COMPARISON_TEST("MEG${ADJOINT}EST-dip${SENSORORIENT}-${HEAD}-dip6"
                        ${HEAD}-dip${SENSORORIENT}.est_meg${ADJOINT} analytic/meg_analytic${SENSORORIENT}.txt -eps 0.001 -col 6 -full
                        DEPENDS MEG${ADJOINT}-dipoles-${HEAD}${SENSORORIENT})
                ENDFOREACH()
            ENDFOREACH()
        ENDIF()
    ENDFOREACH()
ENDFOREACH()

# set tests that are expected to fail :
FOREACH(HEADGEO ${NNc1})
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip4-tangential-rdm PROPERTIES WILL_FAIL TRUE)
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
ENDFOREACH()
FOREACH(HEADGEO "1" ${NNa1} ${NNb1})
    FOREACH(DIP 2 3 4 5)
        FOREACH(COMP mag rdm)
            SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip${DIP}-tangential-${COMP} PROPERTIES WILL_FAIL TRUE)
        ENDFOREACH()
    ENDFOREACH()
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip4-noradial-mag   PROPERTIES WILL_FAIL TRUE)
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip5-noradial-mag   PROPERTIES WILL_FAIL TRUE)
ENDFOREACH()
FOREACH(HEADGEO ${NNc2})
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
ENDFOREACH()
FOREACH(HEADGEO "2" ${NNa2} ${NNb2})
    FOREACH(DIP 4 5)
        FOREACH(COMP mag rdm)
            IF (NOT(${HEADGEO} STREQUAL "NNa2") OR NOT(${DIP} STREQUAL 4) OR NOT(${COMP} STREQUAL mag))
                SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip${DIP}-tangential-${COMP} PROPERTIES WILL_FAIL TRUE)
            ENDIF()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()
