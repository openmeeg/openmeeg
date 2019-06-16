include(ComparisonTest.cmake)

set_file_properties(Suffix
    "HM .hm" "HMInv .hm_inv" "DSM .dsm" "SSM .ssm" "H2EM .h2em" "SurfGainEEG .sgem"
    "ESTEEG .est_eeg" "EEGadjointEST .est_eegadjoint" "MEGadjointEST .est_megadjoint"
    "H2MM .h2mm" "SS2MM .ss2mm" "SurfGainMEG .sgmm" "ESTMEG .est_meg"
)

set_file_properties(CompareOptions
    "HM -sym" "HMInv -sym" "DSM -full:-if1:binary:-if2:binary" "SSM -full:-if1:binary:-if2:binary" "H2EM -if:binary:-sparse"
    "SurfGainEEG -full" "ESTEEG -full" "H2MM -full" "SS2MM -full" "SurfGainMEG -full" "ESTMEG -full"
)


#   TEST COMMON RESULTS ON HEAD1 (Regression test)

foreach (COMPARISON HM HMInv SSM DSM H2EM SurfGainEEG ESTEEG H2MM SS2MM SurfGainMEG ESTMEG)
    set(BASE_FILE_NAME Head1${Suffix_${COMPARISON}})
    OPENMEEG_COMPARISON_TEST(${COMPARISON}-Head1 ${BASE_FILE_NAME} initialTest/${BASE_FILE_NAME} ${CompareOptions_${COMPARISON}})
endforeach()

# Verify ECoG transfert matrices.
# Verify that old and new call for H2ECOGM provide the same answer.

foreach (HEADNUM 1 2 ${HEAD3})
    set(BASE_FILE_NAME Head${HEADNUM}.ecog)
    OPENMEEG_COMPARISON_TEST(H2ECOGM-Head${HEADNUM} ${BASE_FILE_NAME} initialTest/${BASE_FILE_NAME} "-sparse")
    OPENMEEG_COMPARISON_TEST(H2ECOGM-OLD-Head${HEADNUM} Head${HEADNUM}-old.ecog initialTest/${BASE_FILE_NAME} "-sparse")
endforeach()

#   TEST EEG RESULTS ON DIPOLES

# defining variables for those who do not use VTK
# TODO handle comparison test for NNX1 geom when not USE_VTK...

if (USE_VTK) 
    set(NNG "NNa" "NNb" "NNc")
    set(NNGa "NNa")
    set(NNGb "NNb")
    set(NNGc "NNc")
    set(NNGa1 "NNa1")
    set(NNGb1 "NNb1")
    set(NNGc1 "NNc1")
    set(NNGa2 "NNa2")
    set(NNGb2 "NNb2")
    set(NNGc2 "NNc2")
endif()

if (TEST_HEAD3)
    set(HEAD3 3)
endif()

set(EPSILON1 0.13)
set(EPSILON2 0.1)
set(EPSILON3 0.03)
foreach(DIP 1 2 3 4 5)
    foreach(HEADGEO "" ${NNG})
        foreach(HEADNUM 1 2 ${HEAD3})
            foreach(COMP mag rdm)
                set(HEAD "Head${HEADGEO}${HEADNUM}")
                foreach(ADJOINT "" adjoint adjoint2)
                    set(BASE_FILE_NAME "${HEAD}-dip.est_eeg${ADJOINT}")
                    # Compare EEG result with analytical solution obtained with Matlab
                    OPENMEEG_COMPARISON_TEST("EEG${ADJOINT}EST-dip-${HEAD}-dip${DIP}-${COMP}"
                        ${BASE_FILE_NAME} analytic/eeg_head${HEADNUM}_analytic.txt -${COMP} -eps ${EPSILON${HEADNUM}} -col ${DIP} -full
                        DEPENDS EEG${ADJOINT}-dipoles-${HEAD})
                endforeach()
            endforeach()
        endforeach()
    endforeach()
endforeach()
set(EPSILON 0.13)
if (TEST_HEAD3)
    foreach(DIP 1 2)
        foreach(COMP mag rdm)
            set(BASE_FILE_NAME "Head3-dip-skullscalp.est_eeg")
            # Compare EEG result with a numerical solution obtained with the tetrahedral FEM
            OPENMEEG_COMPARISON_TEST("EEGEST-dipSkullScalp-Head3-dip${DIP}-${COMP}"
                ${BASE_FILE_NAME} analytic/eeg_head3-skullscalp_FEM.txt -${COMP} -eps ${EPSILON} -col ${DIP} -full
                DEPENDS EEG-dipolesSkullScalp-Head3)
        endforeach()
    endforeach()
endif()

#   Set tests that are expected to fail :

foreach(ADJOINT "" "adjoint" "adjoint2")
    foreach(HEADGEO ${NNGc1})
        foreach(DIP 1 2 3 4 5)
            set_tests_properties(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip${DIP}-mag PROPERTIES WILL_FAIL TRUE) # all cmp-EEG-mag NNGc1 tests fail...
        endforeach()
    endforeach()
    foreach(HEADGEO "1" ${NNGa1} ${NNGb1})
        set_tests_properties(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip4-rdm PROPERTIES WILL_FAIL TRUE)
        set_tests_properties(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip5-rdm PROPERTIES WILL_FAIL TRUE)
    endforeach()
endforeach()

set(VALIDATION_EIT "${CMAKE_CURRENT_BINARY_DIR}/../tests/test_validationEIT")
if (WIN32)
    set(VALIDATION_EIT "${EXECUTABLE_OUTPUT_PATH}/../tests/test_validationEIT")
endif()

# EIT tests
foreach(HEADGEO "")
    foreach(HEADNUM 1 2 ${HEAD3})
        set(HEAD "Head${HEADGEO}${HEADNUM}")
        set(DATADIR ${OpenMEEG_SOURCE_DIR}/data)
        set(MODELPREFIX ${DATADIR}/${HEAD}/${HEAD})
        set(TESTPREFIX ${OpenMEEG_BINARY_DIR}/tests/${HEAD})
        # validationEIT geometry.geom conductivity.cond dipoles.dip source.dsm headmatinv.bin out.eit_qgradVj out.eit_diffVf
        OPENMEEG_TEST(
            EITvalidation-dipoles-${HEAD}
            ${VALIDATION_EIT} ${MODELPREFIX}.geom ${MODELPREFIX}.cond
            ${DATADIR}/${HEAD}/${HEAD}.dip
            ${TESTPREFIX}.dsm ${TESTPREFIX}.hm_inv ${TESTPREFIX}.eit_qgradVj ${TESTPREFIX}.eit_diffVf
            DEPENDS HMInv-${HEAD} DSM-${HEAD})
        foreach(DIP 1 2 3 4 5)
             # Compare the q.gradVj to diffVf in order to validate the EIT problem
            OPENMEEG_COMPARISON_TEST("EIT-${HEAD}-dip${DIP}"
                ${HEAD}.eit_qgradVj ${OpenMEEG_BINARY_DIR}/tests/${HEAD}.eit_diffVf -eps ${EPSILON${HEADNUM}} -col ${DIP} -full
                DEPENDS EITvalidation-dipoles-${HEAD})
        endforeach()
    endforeach()
endforeach()

foreach(DIP 1 2 3 4 5)
    foreach(HEADGEO "" ${NNG})
        foreach(HEADNUM 1 2 ${HEAD3})
            set(HEAD "Head${HEADGEO}${HEADNUM}")
            # Compare the potential results in a interior sphere of the Surf2Vol operator with analytical solution
            # obtained with Sphere (V.Hedou Modified)
             OPENMEEG_COMPARISON_TEST("EEGinternal-dip-${HEAD}-dip${DIP}"
                                      ${HEAD}-dip-internal.est_eeg analytic/eeg_internal_analytic.txt -eps ${EPSILON${HEADNUM}} -col ${DIP} -full
                                      DEPENDS InternalPot-dipoles-${HEAD})
        endforeach()
    endforeach()
endforeach()

#   ECoG comparison tests
#   Done only on dipoles 1, 2 and 3 because otherwise the error (even in the analytic model is poorly controlled.

set(DIP_EPSILON1 0.035)
set(DIP_EPSILON2 0.17)
set(DIP_EPSILON3 0.59)
foreach(HEADNUM 1 2 ${HEAD3})
    foreach(DIP 1 2 3)
        OPENMEEG_COMPARISON_TEST("GAINECOG-Head${HEADNUM}-dip${DIP}"
                                 Head${HEADNUM}ECoGGain.mat analytic/ecog_analytic.txt -eps ${DIP_EPSILON${DIP}} -col ${DIP} -full
                                 DEPENDS GAINECOG-Head${HEADNUM})
    endforeach()
endforeach()

#   Set tests that are expected to fail :

foreach(HEADGEO "1" ${NNGa1} ${NNGb1})
    set_tests_properties(cmp-EEGinternal-dip-Head${HEADGEO}-dip5 PROPERTIES WILL_FAIL TRUE)
endforeach()

set(EPSILON1 0.15)
set(EPSILON2 0.14)
set(EPSILON3 0.09)

foreach(SENSORORIENT "" "-tangential" "-noradial")
    foreach(ADJOINT "" adjoint adjoint2)
        if (NOT(${ADJOINT} STREQUAL "adjoint" OR ${ADJOINT} STREQUAL "adjoint2") OR (SENSORORIENT STREQUAL ""))
            foreach(HEADGEO "" ${NNG})
                foreach(HEADNUM 1 2 ${HEAD3})
                    set(HEAD "Head${HEADGEO}${HEADNUM}")
                    foreach(DIP 1 2 3 4 5)
                        foreach(COMP mag rdm)
                            # Compare MEG result with analytical solution obtained with Matlab
                            OPENMEEG_COMPARISON_TEST("MEG${ADJOINT}EST-dip-${HEAD}-dip${DIP}${SENSORORIENT}-${COMP}"
                                ${HEAD}-dip${SENSORORIENT}.est_meg${ADJOINT} analytic/meg_analytic${SENSORORIENT}.txt -${COMP} -eps ${EPSILON${HEADNUM}} -col ${DIP} -full
                                DEPENDS MEG${ADJOINT}-dipoles-${HEAD}${SENSORORIENT})
                        endforeach()
                    endforeach()
                    OPENMEEG_COMPARISON_TEST("MEG${ADJOINT}EST-dip${SENSORORIENT}-${HEAD}-dip6"
                        ${HEAD}-dip${SENSORORIENT}.est_meg${ADJOINT} analytic/meg_analytic${SENSORORIENT}.txt -eps 0.001 -col 6 -full
                        DEPENDS MEG${ADJOINT}-dipoles-${HEAD}${SENSORORIENT})
                endforeach()
            endforeach()
        endif()
    endforeach()
endforeach()

foreach(HEADGEO ${NNGc1})
    set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip4-tangential-rdm PROPERTIES WILL_FAIL TRUE)
    set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(HEADGEO "1" ${NNGa1} ${NNGb1})
    foreach(DIP 2 3 4 5)
        foreach(COMP mag rdm)
            set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip${DIP}-tangential-${COMP} PROPERTIES WILL_FAIL TRUE)
        endforeach()
    endforeach()
    set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip4-noradial-mag   PROPERTIES WILL_FAIL TRUE)
    set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip5-noradial-mag   PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(HEADGEO ${NNGc2})
    set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(HEADGEO "2" ${NNGa2} ${NNGb2})
    foreach(DIP 4 5)
        foreach(COMP mag rdm)
            if (NOT(${HEADGEO} STREQUAL "NNa2") OR NOT(${DIP} STREQUAL 4) OR NOT(${COMP} STREQUAL mag))
                set_tests_properties(cmp-MEGEST-dip-Head${HEADGEO}-dip${DIP}-tangential-${COMP} PROPERTIES WILL_FAIL TRUE)
            endif()
        endforeach()
    endforeach()
endforeach()

# tests for the multiple nonconductive branch
# compare with analytical solution
set(EPSILONMN1 0.10)
set(EPSILONMN2 0.10)
foreach(DIP 1 2 3 4 5 6)
    foreach(HEADNUM MN1 MN2)
        foreach(COMP mag rdm)
        set(HEAD "Head${HEADNUM}")
            foreach(ADJOINT "" adjoint adjoint2)
                set(BASE_FILE_NAME "${HEAD}-dip.est_eeg${ADJOINT}")
                OPENMEEG_COMPARISON_TEST("EEG${ADJOINT}EST-dip-${HEAD}-dip${DIP}-${COMP}" ${BASE_FILE_NAME} analytic/eeg_head${HEADNUM}_analytic.txt 
                    -${COMP} -eps ${EPSILON${HEADNUM}} -col ${DIP} -full DEPENDS EEG${ADJOINT}-dipoles-${HEAD})
            endforeach()
        endforeach()
    endforeach()
endforeach()
