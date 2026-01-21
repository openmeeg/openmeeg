#   This is the only test that should not use OPENMEEG_TEST

add_test(
    NAME CLEAN-TESTS
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/delete_head_files_from_working_dir.cmake
)

set(APPS_DIR ${PROJECT_BINARY_DIR}/apps)
set(ASSEMBLE ${APPS_DIR}/om_assemble)
set(INVERSER ${APPS_DIR}/om_minverser)
set(GAIN     ${APPS_DIR}/om_gain)
set(FORWARD  ${APPS_DIR}/om_forward)

OPENMEEG_TEST(assemble-help ${ASSEMBLE} -h)
OPENMEEG_TEST(inverser-help ${INVERSER} -h)
OPENMEEG_TEST(gain-help ${GAIN} -h)
OPENMEEG_TEST(forward-help ${FORWARD} -h)

function(TESTHEAD HEADNUM)
    set(SUBJECT "Head${HEADNUM}")
    set(MODELBASE ${PROJECT_SOURCE_DIR}/data/${SUBJECT}/${SUBJECT})

    set(GEOM                   ${MODELBASE}.geom)
    set(COND                   ${MODELBASE}.cond)
    set(SRCMESH                ${MODELBASE}.tri)
    set(POINTS                 ${CMAKE_CURRENT_SOURCE_DIR}/analytic/eeg_internal_points.txt)
    set(MONOPPOS               ${MODELBASE}.monop)
    set(DIPPOS                 ${MODELBASE}.dip)
    set(DIPPOS-SKULLSCALP      ${MODELBASE}-skullscalp.dip)
    set(PATCHES                ${MODELBASE}.patches)
    set(EITPATCHES             ${MODELBASE}-EIT.patches)
    set(SQUIDS                 ${MODELBASE}.squids)
    set(ECOG-ELECTRODES        ${MODELBASE}-ecog.electrodes)
    set(SQUIDS-TANGENTIAL      ${MODELBASE}-tangential.squids)
    set(SQUIDS-NORADIAL        ${MODELBASE}-noradial.squids)
    set(SURFSOURCES            ${MODELBASE}.src)
    set(DIPSOURCES             ${MODELBASE}.srcdip)
    set(DIPSOURCES-SKULLSCALP  ${MODELBASE}-skullscalp.srcdip)

    set(SMOOTH                 ${SUBJECT}.smooth)
    set(AREAS                  ${SUBJECT}.ai)
    set(HMMAT                  ${SUBJECT}.hm)
    set(HMINVMAT               ${SUBJECT}.hm_inv)
    set(SSMMAT                 ${SUBJECT}.ssm)
    set(CMMAT                  ${SUBJECT}.cm)
    set(ECOGMMAT               ${SUBJECT}.ecog)
    set(OLDECOGMMAT            ${SUBJECT}-old.ecog)
    set(GAINECOGMAT            ${SUBJECT}ECoGGain.mat)
    set(H2EMMAT                ${SUBJECT}.h2em)
    set(SGEMMAT                ${SUBJECT}.sgem)
    set(H2IPMAT                ${SUBJECT}.h2ip)
    set(EITSMAT                ${SUBJECT}.eitsm)
    set(H2MMMAT                ${SUBJECT}.h2mm)
    set(H2MMMAT-TANGENTIAL     ${SUBJECT}-tangential.h2mm)
    set(H2MMMAT-NORADIAL       ${SUBJECT}-noradial.h2mm)
    set(SS2MMMAT               ${SUBJECT}.ss2mm)
    set(SGMMMAT                ${SUBJECT}.sgmm)
    set(DS2IPMAT               ${SUBJECT}.ds2ip)
    set(MSMMAT                 ${SUBJECT}.msm)
    set(DSMMAT                 ${SUBJECT}.dsm)
    set(DSM-SKULLSCALPMAT      ${SUBJECT}-skullscalp.dsm)
    set(DS2MMMAT               ${SUBJECT}.ds2mm)
    set(DS2MMMAT-TANGENTIAL    ${SUBJECT}-tangential.ds2mm)
    set(DS2MMMAT-NORADIAL      ${SUBJECT}-noradial.ds2mm)
    set(DGIPMAT                ${SUBJECT}.dgip)
    set(GSIPMAT                ${SUBJECT}.gsip)
    set(DGEMMAT                ${SUBJECT}.dgem)
    set(DGEM-SKULLSCALPMAT     ${SUBJECT}-skullscalp.dgem)
    set(DGEMADJOINTMAT         ${SUBJECT}-adjoint.dgem)
    set(DGEMADJOINT2MAT        ${SUBJECT}-adjoint2.dgem)
    set(DGMMMAT                ${SUBJECT}.dgmm)
    set(DGMMADJOINTMAT         ${SUBJECT}-adjoint.dgmm)
    set(DGMMADJOINT2MAT        ${SUBJECT}-adjoint2.dgmm)
    set(DGMMMAT-TANGENTIAL     ${SUBJECT}-tangential.dgmm)
    set(DGMMMAT-NORADIAL       ${SUBJECT}-noradial.dgmm)

    set(ESTEEG                 ${SUBJECT}.est_eeg)
    set(ESTMEG                 ${SUBJECT}.est_meg)
    set(EEGINVBASE             ${SUBJECT}-eeg)
    set(MEGINVBASE             ${SUBJECT}-meg)
    set(ESTDIPBASE             ${SUBJECT}-dip)

    #   assemble tests

    OPENMEEG_TEST(HM-${SUBJECT} ${ASSEMBLE} -HM ${GEOM} ${COND} ${HMMAT} DEPENDS CLEAN-TESTS)
    OPENMEEG_TEST(HMInv-${SUBJECT} ${INVERSER} ${HMMAT} ${HMINVMAT}      DEPENDS HM-${SUBJECT})

    if (${HEADNUM} EQUAL 1)

        OPENMEEG_TEST(SSM-${SUBJECT} ${ASSEMBLE} -SSM ${GEOM} ${COND} ${SRCMESH} ${SSMMAT} DEPENDS CLEAN-TESTS)

        # corticalMat tests

        OPENMEEG_TEST(CM1-0-${SUBJECT} ${ASSEMBLE} -CM ${GEOM} ${COND} ${PATCHES} "Brain" ${CMMAT} DEPENDS CLEAN-TESTS)
        OPENMEEG_TEST(CM1-1-${SUBJECT} ${ASSEMBLE} -CM ${GEOM} ${COND} ${PATCHES} "Brain" ${CMMAT} 1e-4 1.58e-2 DEPENDS CLEAN-TESTS)
        OPENMEEG_TEST(CM2-${SUBJECT}   ${ASSEMBLE} -CM ${GEOM} ${COND} ${PATCHES} "Brain" ${CMMAT} 12.4 DEPENDS CLEAN-TESTS)

    endif()

    # EEG test

    OPENMEEG_TEST(H2EM-${SUBJECT} ${ASSEMBLE} -H2EM ${GEOM} ${COND} ${PATCHES} ${H2EMMAT} DEPENDS CLEAN-TESTS)

    if (${HEADNUM} EQUAL 1)

        OPENMEEG_TEST(SurfGainEEG-${SUBJECT} ${GAIN} -EEG ${HMINVMAT} ${SSMMAT} ${H2EMMAT} ${SGEMMAT}
                      DEPENDS HMInv-${SUBJECT} SSM-${SUBJECT} H2EM-${SUBJECT})

        # Forward

        OPENMEEG_TEST(ESTEEG-${SUBJECT} ${FORWARD} ${SGEMMAT} ${SURFSOURCES} ${ESTEEG} 0.0 DEPENDS SurfGainEEG-${SUBJECT})

    endif()

    # MEG

    OPENMEEG_TEST(H2MM-${SUBJECT} ${ASSEMBLE} -H2MM ${GEOM} ${COND} ${SQUIDS} ${H2MMMAT} DEPENDS CLEAN-TESTS)
    OPENMEEG_TEST(H2MM-${SUBJECT}-tangential ${ASSEMBLE} -H2MM ${GEOM} ${COND} ${SQUIDS-TANGENTIAL} ${H2MMMAT-TANGENTIAL} DEPENDS CLEAN-TESTS)
    OPENMEEG_TEST(H2MM-${SUBJECT}-noradial ${ASSEMBLE} -H2MM ${GEOM} ${COND} ${SQUIDS-NORADIAL} ${H2MMMAT-NORADIAL} DEPENDS CLEAN-TESTS)

    if (${HEADNUM} EQUAL 1)

        OPENMEEG_TEST(SS2MM-${SUBJECT} ${ASSEMBLE} -SS2MM ${SRCMESH} ${SQUIDS} ${SS2MMMAT} DEPENDS CLEAN-TESTS)

        OPENMEEG_TEST(SurfGainMEG-${SUBJECT} ${GAIN} -MEG ${HMINVMAT} ${SSMMAT} ${H2MMMAT} ${SS2MMMAT} ${SGMMMAT}
                      DEPENDS HMInv-${SUBJECT} SSM-${SUBJECT} H2MM-${SUBJECT} SS2MM-${SUBJECT})

        # Forward
        OPENMEEG_TEST(ESTMEG-${SUBJECT} ${FORWARD} ${SGMMMAT} ${SURFSOURCES} ${ESTMEG} 0.0
                      DEPENDS SurfGainMEG-${SUBJECT})
    endif()

    # ECoG

    if (${HEADNUM} EQUAL 1 OR ${HEADNUM} EQUAL 2 OR ${HEADNUM} EQUAL 3)
        OPENMEEG_TEST(H2ECOGM-OLD-${SUBJECT} ${ASSEMBLE} -H2ECOGM ${GEOM} ${COND} ${ECOG-ELECTRODES} ${OLDECOGMMAT} DEPENDS CLEAN-TESTS)
        OPENMEEG_TEST(H2ECOGM-${SUBJECT} ${ASSEMBLE} -H2ECOGM ${GEOM} ${COND} ${ECOG-ELECTRODES} "1" ${ECOGMMAT} DEPENDS CLEAN-TESTS)
    endif()

    # TEST DIPOLE FORWARD RESULTS (Regression test)

    # om_assemble -MSM geometry.geom conductivity.cond monooples.monop msm.bin

    OPENMEEG_TEST(MSM-${SUBJECT} ${ASSEMBLE} -MSM ${GEOM} ${COND} ${MONOPPOS} ${MSMMAT} DEPENDS CLEAN-TESTS)

    # om_assemble -DSM geometry.geom conductivity.cond dipoles.dip dsm.bin

    OPENMEEG_TEST(DSM-${SUBJECT} ${ASSEMBLE} -DSM ${GEOM} ${COND} ${DIPPOS} ${DSMMAT} DEPENDS CLEAN-TESTS)

    # om_assemble -DS2MM dipoles.dip squidscoord.squids sToMEGmat.bin

    OPENMEEG_TEST(DS2MM-${SUBJECT} ${ASSEMBLE} -DS2MM ${DIPPOS} ${SQUIDS} ${DS2MMMAT} DEPENDS CLEAN-TESTS)
    OPENMEEG_TEST(DS2MM-${SUBJECT}-tangential ${ASSEMBLE} -DS2MM ${DIPPOS} ${SQUIDS-TANGENTIAL} ${DS2MMMAT-TANGENTIAL} DEPENDS CLEAN-TESTS)
    OPENMEEG_TEST(DS2MM-${SUBJECT}-noradial ${ASSEMBLE} -DS2MM ${DIPPOS} ${SQUIDS-NORADIAL} ${DS2MMMAT-NORADIAL} DEPENDS CLEAN-TESTS)

    #   ECoG gain matrix.

    if (${HEADNUM} EQUAL 1 OR ${HEADNUM} EQUAL 2 OR ${HEADNUM} EQUAL 3)
        OPENMEEG_TEST(GAINECOG-${SUBJECT} ${GAIN} -EEG ${HMINVMAT} ${DSMMAT} ${ECOGMMAT} ${GAINECOGMAT} DEPENDS HMInv-${SUBJECT} DSM-${SUBJECT} H2ECOGM-${SUBJECT})
    endif()

    # om_assemble -H2IPM geometry.geom conductivity.cond internalpoints.txt h2ip.bin

    OPENMEEG_TEST(H2IPM-${SUBJECT} ${ASSEMBLE} -H2IPM ${GEOM} ${COND} ${POINTS} ${H2IPMAT} DEPENDS CLEAN-TESTS)

    # for Head1 and Head2, test EIT
    if (${HEADNUM} EQUAL 1 OR ${HEADNUM} EQUAL 2)
        # EIT InternalPot
        OPENMEEG_TEST(EITSM-${SUBJECT} ${ASSEMBLE} -EITSM ${GEOM} ${COND} ${EITPATCHES} ${EITSMAT})

        OPENMEEG_TEST(GainEITInternalPot-${SUBJECT} ${GAIN} -EITIP ${HMINVMAT} ${EITSMAT} ${H2IPMAT} ${GSIPMAT}
                      DEPENDS HMInv-${SUBJECT} EITSM-${SUBJECT} H2IPM-${SUBJECT})
    endif()

    # om_assemble -DS2IPM geometry.geom conductivity.cond dipoles.dip internalpoints.txt ds2ip.bin

    OPENMEEG_TEST(S2IPM-${SUBJECT} ${ASSEMBLE} -DS2IPM ${GEOM} ${COND} ${DIPPOS} ${POINTS} ${DS2IPMAT} DEPENDS CLEAN-TESTS)

    OPENMEEG_TEST(DipGainEEG-${SUBJECT} ${GAIN} -EEG ${HMINVMAT} ${DSMMAT} ${H2EMMAT} ${DGEMMAT}
                  DEPENDS HMInv-${SUBJECT} DSM-${SUBJECT} H2EM-${SUBJECT})
    OPENMEEG_TEST(DipGainEEGadjoint-${SUBJECT} ${GAIN} -EEGadjoint ${GEOM} ${COND} ${DIPPOS} ${HMMAT} ${H2EMMAT} ${DGEMADJOINTMAT}
                  DEPENDS HM-${SUBJECT} H2EM-${SUBJECT})
    OPENMEEG_TEST(DipGainMEG-${SUBJECT} ${GAIN} -MEG ${HMINVMAT} ${DSMMAT} ${H2MMMAT} ${DS2MMMAT} ${DGMMMAT}
                  DEPENDS HMInv-${SUBJECT} DSM-${SUBJECT} H2MM-${SUBJECT} DS2MM-${SUBJECT})
    OPENMEEG_TEST(DipGainMEGadjoint-${SUBJECT} ${GAIN} -MEGadjoint ${GEOM} ${COND} ${DIPPOS} ${HMMAT} ${H2MMMAT} ${DS2MMMAT} ${DGMMADJOINTMAT}
                  DEPENDS HM-${SUBJECT} H2MM-${SUBJECT} DS2MM-${SUBJECT})
    OPENMEEG_TEST(DipGainMEG-${SUBJECT}-tangential ${GAIN} -MEG ${HMINVMAT} ${DSMMAT} ${H2MMMAT-TANGENTIAL} ${DS2MMMAT-TANGENTIAL} ${DGMMMAT-TANGENTIAL}
                  DEPENDS HMInv-${SUBJECT} DSM-${SUBJECT} H2MM-${SUBJECT}-tangential DS2MM-${SUBJECT}-tangential)
    OPENMEEG_TEST(DipGainMEG-${SUBJECT}-noradial ${GAIN} -MEG ${HMINVMAT} ${DSMMAT} ${H2MMMAT-NORADIAL} ${DS2MMMAT-NORADIAL} ${DGMMMAT-NORADIAL}
                  DEPENDS HMInv-${SUBJECT} DSM-${SUBJECT} H2MM-${SUBJECT}-noradial DS2MM-${SUBJECT}-noradial)
    OPENMEEG_TEST(DipGainEEGMEGadjoint-${SUBJECT} ${GAIN} -EEGMEGadjoint ${GEOM} ${COND} ${DIPPOS} ${HMMAT} ${H2EMMAT} ${H2MMMAT} ${DS2MMMAT} ${DGEMADJOINT2MAT} ${DGMMADJOINT2MAT}
                  DEPENDS HM-${SUBJECT} H2EM-${SUBJECT} H2MM-${SUBJECT} DS2MM-${SUBJECT})
    OPENMEEG_TEST(DipGainInternalPot-${SUBJECT} ${GAIN} -IP ${HMINVMAT} ${DSMMAT} ${H2IPMAT} ${DS2IPMAT} ${DGIPMAT}
                  DEPENDS HMInv-${SUBJECT} DSM-${SUBJECT} H2IPM-${SUBJECT} S2IPM-${SUBJECT})

    # forward gainmatrix.bin dipoleActivation.src estimatedeegdata.txt noiselevel

    OPENMEEG_TEST(EEG-dipoles-${SUBJECT} ${FORWARD} ${DGEMMAT} ${DIPSOURCES} ${ESTDIPBASE}.est_eeg 0.0
                  DEPENDS DipGainEEG-${SUBJECT})
    OPENMEEG_TEST(EEGadjoint-dipoles-${SUBJECT} ${FORWARD} ${DGEMADJOINTMAT} ${DIPSOURCES} ${ESTDIPBASE}.est_eegadjoint 0.0
                  DEPENDS DipGainEEGadjoint-${SUBJECT})
    OPENMEEG_TEST(EEGadjoint2-dipoles-${SUBJECT} ${FORWARD} ${DGEMADJOINT2MAT} ${DIPSOURCES} ${ESTDIPBASE}.est_eegadjoint2 0.0
                  DEPENDS DipGainEEGMEGadjoint-${SUBJECT})
    OPENMEEG_TEST(MEG-dipoles-${SUBJECT} ${FORWARD} ${DGMMMAT} ${DIPSOURCES} ${ESTDIPBASE}.est_meg 0.0
                  DEPENDS DipGainMEG-${SUBJECT})
    OPENMEEG_TEST(MEGadjoint-dipoles-${SUBJECT} ${FORWARD} ${DGMMADJOINTMAT} ${DIPSOURCES} ${ESTDIPBASE}.est_megadjoint 0.0
                  DEPENDS DipGainMEGadjoint-${SUBJECT})
    OPENMEEG_TEST(MEG-dipoles-${SUBJECT}-tangential ${FORWARD} ${DGMMMAT-TANGENTIAL} ${DIPSOURCES} ${ESTDIPBASE}-tangential.est_meg 0.0
                  DEPENDS DipGainMEG-${SUBJECT}-tangential)
    OPENMEEG_TEST(MEG-dipoles-${SUBJECT}-noradial ${FORWARD} ${DGMMMAT-NORADIAL} ${DIPSOURCES} ${ESTDIPBASE}-noradial.est_meg 0.0
                  DEPENDS DipGainMEG-${SUBJECT}-noradial)
    OPENMEEG_TEST(MEGadjoint2-dipoles-${SUBJECT} ${FORWARD} ${DGMMADJOINT2MAT} ${DIPSOURCES} ${ESTDIPBASE}.est_megadjoint2 0.0
                  DEPENDS DipGainEEGMEGadjoint-${SUBJECT})
    OPENMEEG_TEST(InternalPot-dipoles-${SUBJECT} ${FORWARD} ${DGIPMAT} ${DIPSOURCES} ${ESTDIPBASE}-internal.est_eeg 0.0
                  DEPENDS DipGainInternalPot-${SUBJECT})

    # tests on Head3 for dipoles in the skull and scalp
    if (${HEADNUM} EQUAL 3)
        OPENMEEG_TEST(DSMSkullScalp-${SUBJECT} ${ASSEMBLE} -DSM ${GEOM} ${COND} ${DIPPOS-SKULLSCALP} ${DSM-SKULLSCALPMAT} DEPENDS CLEAN-TESTS)

        OPENMEEG_TEST(DipGainEEGSkullScalp-${SUBJECT} ${GAIN} -EEG ${HMINVMAT} ${DSM-SKULLSCALPMAT} ${H2EMMAT} ${DGEM-SKULLSCALPMAT}
                DEPENDS HMInv-${SUBJECT} DSMSkullScalp-${SUBJECT} H2EM-${SUBJECT})

        OPENMEEG_TEST(EEG-dipolesSkullScalp-${SUBJECT} ${FORWARD} ${DGEM-SKULLSCALPMAT} ${DIPSOURCES-SKULLSCALP} ${ESTDIPBASE}-skullscalp.est_eeg 0.0
                DEPENDS DipGainEEGSkullScalp-${SUBJECT})
    endif()
endfunction()
