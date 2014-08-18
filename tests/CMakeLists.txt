include(ComparisonTest)

set_file_properties(Suffix
    "HM .hm" "HMInv .hm_inv" "DSM .dsm" "SSM .ssm" "v2eeg .h2em" "SurfGainEEG .sgem"
    "ESTEEG .est_eeg" "EEGadjointEST .est_eegadjoint" "MEGadjointEST .est_megadjoint"
    "Head2MEGMat .h2mm" "SS2MM .ss2mm" "SurfGainMEG .sgmm" "ESTMEG .est_meg"
)

set_file_properties(CompareOptions
    "HM -sym" "HMInv -sym" "DSM -full:-if1:binary:-if2:binary" "SSM -full:-if1:binary:-if2:binary" "v2eeg -if:binary:-sparse"
    "SurfGainEEG -full" "ESTEEG -full" "Head2MEGMat -full" "SS2MM -full" "SurfGainMEG -full" "ESTMEG -full"
)

#   IO TESTS

OPENMEEG_UNIT_TEST(load_geo
    SOURCES load_geo.cpp
    LIBRARIES OpenMEEG OpenMEEGMaths
    PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.cond)

OPENMEEG_UNIT_TEST(test_mesh_ios
    SOURCES test_mesh_ios.cpp
    LIBRARIES OpenMEEG OpenMEEGMaths
    PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.tri)

NEW_EXECUTABLE(test_sensors test_sensors.cpp
               LIBRARIES OpenMEEG)
NEW_EXECUTABLE(compare_matrix compare_matrix.cpp
               LIBRARIES OpenMEEGMaths ${LAPACK_LIBRARIES})
NEW_EXECUTABLE(om_validationEIT validationEIT.cpp ${OPEMEEG_HEADERS}
               LIBRARIES OpenMEEG OpenMEEGMaths ${LAPACK_LIBRARIES})

if (WIN32 OR APPLE)
    set(INSTALLED_TESTS ${INSTALLED_TESTS} test_sensors compare_matrix om_validationEIT)
    install(TARGETS test_sensors compare_matrix om_validationEIT
            RUNTIME DESTINATION ${INSTALL_BIN_DIR}
            LIBRARY DESTINATION ${INSTALL_LIB_DIR}
            ARCHIVE DESTINATION ${INSTALL_LIB_DIR})
endif()

#   TEST COMMON RESULTS ON HEAD1 (Regression test)

set(COMPARISONS HM HMInv SSM DSM v2eeg SurfGainEEG ESTEEG
                Head2MEGMat SS2MM SurfGainMEG ESTMEG)

foreach (COMPARISON ${COMPARISONS})
    set(BASE_FILE_NAME Head1${Suffix_${COMPARISON}})
    set(TEST_NAME ${COMPARISON}-Head1)
    OPENMEEG_COMPARISON_TEST(${COMPARISON}-Head1
        ${BASE_FILE_NAME} initialTest/${BASE_FILE_NAME} ${CompareOptions_${COMPARISON}})
endforeach()

#   TEST EEG RESULTS ON DIPOLES

# defining variables for those who do not use VTK
# TODO handle comparison test for NNX1 geom when not USE_VTK...

if (USE_VTK) 
    set(NN "NNa" "NNb" "NNc")
    set(NNa "NNa")
    set(NNb "NNb")
    set(NNc "NNc")
    set(NNa1 "NNa1")
    set(NNb1 "NNb1")
    set(NNc1 "NNc1")
    set(NNa2 "NNa2")
    set(NNb2 "NNb2")
    set(NNc2 "NNc2")
endif()

if (TEST_HEAD3)
    set(HEAD3 3)
endif()

set(EPSILON1 0.13)
set(EPSILON2 0.1)
set(EPSILON3 0.03)
foreach(DIP 1 2 3 4 5)
    foreach(HEADGEO "" ${NN})
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
    foreach(HEADGEO ${NNc1})
        foreach(DIP 1 2 3 4 5)
            SET_TESTS_PROPERTIES(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip${DIP}-mag PROPERTIES WILL_FAIL TRUE) # all cmp-EEG-mag NNc1 tests fail...
        endforeach()
    endforeach()
    foreach(HEADGEO "1" ${NNa1} ${NNb1})
        SET_TESTS_PROPERTIES(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip4-rdm PROPERTIES WILL_FAIL TRUE)
        SET_TESTS_PROPERTIES(cmp-EEG${ADJOINT}EST-dip-Head${HEADGEO}-dip5-rdm PROPERTIES WILL_FAIL TRUE)
    endforeach()
endforeach()

set(VALIDATION_EIT "om_validationEIT")

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
            DEPENDS HMINV-${HEAD} DSM-${HEAD})
        foreach(DIP 1 2 3 4 5)
             # Compare the q.gradVj to diffVf in order to validate the EIT problem
            OPENMEEG_COMPARISON_TEST("EIT-${HEAD}-dip${DIP}"
                ${HEAD}.eit_qgradVj ${OpenMEEG_BINARY_DIR}/tests/${HEAD}.eit_diffVf -eps ${EPSILON${HEADNUM}} -col ${DIP} -full
                DEPENDS EITvalidation-dipoles-${HEAD})
        endforeach()
    endforeach()
endforeach()

foreach(DIP 1 2 3 4 5)
    foreach(HEADGEO "" ${NN})
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

#   Set tests that are expected to fail :

foreach(HEADGEO "1" ${NNa1} ${NNb1})
    SET_TESTS_PROPERTIES(cmp-EEGinternal-dip-Head${HEADGEO}-dip5 PROPERTIES WILL_FAIL TRUE)
endforeach()

set(EPSILON1 0.15)
set(EPSILON2 0.14)
set(EPSILON3 0.09)

foreach(SENSORORIENT "" "-tangential" "-noradial")
    foreach(ADJOINT "" adjoint adjoint2)
        if (NOT(${ADJOINT} STREQUAL "adjoint" OR ${ADJOINT} STREQUAL "adjoint2") OR (SENSORORIENT STREQUAL ""))
            foreach(HEADGEO "" ${NN})
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

foreach(HEADGEO ${NNc1})
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip4-tangential-rdm PROPERTIES WILL_FAIL TRUE)
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(HEADGEO "1" ${NNa1} ${NNb1})
    foreach(DIP 2 3 4 5)
        foreach(COMP mag rdm)
            SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip${DIP}-tangential-${COMP} PROPERTIES WILL_FAIL TRUE)
        endforeach()
    endforeach()
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip4-noradial-mag   PROPERTIES WILL_FAIL TRUE)
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip5-noradial-mag   PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(HEADGEO ${NNc2})
    SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip5-tangential-rdm PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(HEADGEO "2" ${NNa2} ${NNb2})
    foreach(DIP 4 5)
        foreach(COMP mag rdm)
            if (NOT(${HEADGEO} STREQUAL "NNa2") OR NOT(${DIP} STREQUAL 4) OR NOT(${COMP} STREQUAL mag))
                SET_TESTS_PROPERTIES(cmp-MEGEST-dip-Head${HEADGEO}-dip${DIP}-tangential-${COMP} PROPERTIES WILL_FAIL TRUE)
            endif()
        endforeach()
    endforeach()
endforeach()
