############ IO TEST ##############
ADD_EXECUTABLE(load_geo load_geo.cpp ${OpenMEEG_SOURCE_DIR}/src/mesh3.cpp)
TARGET_LINK_LIBRARIES (load_geo ${OPENMEEG_OTHER_LIBRARIES})
ADD_TEST(load_tri ${EXECUTABLE_OUTPUT_PATH}/load_geo ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri)

############ ASSEMBLE TEST ##############
ADD_TEST(assemble-LHS ${EXECUTABLE_OUTPUT_PATH}/assemble -LHS ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin)

ADD_TEST(minverser-lhsmatrix ${EXECUTABLE_OUTPUT_PATH}/minverser ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin)

ADD_TEST(assemble-RHS ${EXECUTABLE_OUTPUT_PATH}/assemble -RHS ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin)

ADD_TEST(presmooth ${EXECUTABLE_OUTPUT_PATH}/presmooth ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin)


############ EEG TEST ##############
ADD_TEST(assemble-x2EEGresp ${EXECUTABLE_OUTPUT_PATH}/assemble -xToEEGresponse ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.patches ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin)

ADD_TEST(gain-EEG ${EXECUTABLE_OUTPUT_PATH}/gain -EEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin)

#------ FORWARD -----#
ADD_TEST(forward-EEG ${EXECUTABLE_OUTPUT_PATH}/forward -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.src ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt 0.0)

#------ INVERSE -----#
ADD_TEST(inverse-EEG-HEAT ${EXECUTABLE_OUTPUT_PATH}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGHEAT.txt 1 1e-4 HEAT 10000 0.0)

ADD_TEST(inverse-EEG-MN ${EXECUTABLE_OUTPUT_PATH}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGMN.txt 1 1e-4 MN 10000 0.0)

ADD_TEST(inverse-EEG-TV ${EXECUTABLE_OUTPUT_PATH}/inverse -EEG ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGTV.txt 1 1e-4 TV 10000 0.0)


############ MEG TEST ##############
ADD_TEST(assemble-x2MEGresp ${EXECUTABLE_OUTPUT_PATH}/assemble -xToMEGresponseContrib ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.cond ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin)

ADD_TEST(assemble-s2MEGresp ${EXECUTABLE_OUTPUT_PATH}/assemble -sToMEGresponseContrib ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.squids ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin)

ADD_TEST(gain-MEG ${EXECUTABLE_OUTPUT_PATH}/gain -MEG ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin)

#------ FORWARD -----#
ADD_TEST(forward-MEG ${EXECUTABLE_OUTPUT_PATH}/forward -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_SOURCE_DIR}/Computations/Head1/Head1.src ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt 0.0)

#------ INVERSE -----#
ADD_TEST(inverse-MEG-HEAT ${EXECUTABLE_OUTPUT_PATH}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGHEAT.txt 1 1e-4 HEAT 10000 0.0)

ADD_TEST(inverse-MEG-MN ${EXECUTABLE_OUTPUT_PATH}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGMN.txt 1 1e-4 MN 10000 0.0)

ADD_TEST(inverse-MEG-TV ${EXECUTABLE_OUTPUT_PATH}/inverse -MEG ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/smoothmatrix.bin ${OpenMEEG_BINARY_DIR}/tests/aivector.bin ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGTV.txt 1 1e-4 TV 10000 0.0)


############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

#------ LHS -----#
ADD_TEST(compareLHS ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -sym
                    ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_lhsmatrix.bin)
#------ LHS INV -----#
ADD_TEST(compareLHSINV ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -sym
                       ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin
                       ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_lhsinvmatrix.bin)
#------ RHS -----#
ADD_TEST(compareRHS ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -std
                    ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_rhsmatrix.bin)

############ TEST EEG RESULTS ##############

#------ X2EEG -----#
ADD_TEST(compareX2EEG ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -std 
                      ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_x2eegmatrix.bin)
#------ EEG GAIN -----#
ADD_TEST(compareEEGGAIN ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -std
                        ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_eeggainmatrix.bin)
#------ EEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareEEGEST ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                        ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedeeg.txt)
#------ EEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareEEGINVHEAT ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                           ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGHEAT.txt
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_eegheat.txt)
#------ EEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareEEGINVMN ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGMN.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_eegmn.txt)
#------ EEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareEEGINVTV ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGTV.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_eegtv.txt)

############ TEST MEG RESULTS ##############

#------ X2MEG -----#
ADD_TEST(compareX2MEG ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -std 
                      ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_x2megmatrix.bin)
#------ S2MEG -----#
ADD_TEST(compareS2MEG ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -std
                      ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_s2megmatrix.bin)
#------ MEG GAIN -----#
ADD_TEST(compareMEGGAIN ${EXECUTABLE_OUTPUT_PATH}/comparematrix -bin -std 
                        ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_meggainmatrix.bin)
#------ MEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareMEGEST ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                        ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedmeg.txt)
#------ MEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareMEGINVHEAT ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                           ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGHEAT.txt
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_megheat.txt)
#------ MEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareMEGINVMN ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGMN.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_megmn.txt)
#------ MEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareMEGINVTV ${EXECUTABLE_OUTPUT_PATH}/comparematrix -txt -std
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGTV.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_megtv.txt)
