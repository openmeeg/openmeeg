
ADD_EXECUTABLE(matlibtest test.cpp chrono.cpp)
TARGET_LINK_LIBRARIES (matlibtest ${OPENMEEG_OTHER_LIBRARIES})
ADD_TEST(matlibtest ${EXECUTABLE_OUTPUT_PATH}/matlibtest)

ADD_EXECUTABLE(comparematrix compareMatrix.cpp)
TARGET_LINK_LIBRARIES (comparematrix ${OPENMEEG_OTHER_LIBRARIES})

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
                         
                         
