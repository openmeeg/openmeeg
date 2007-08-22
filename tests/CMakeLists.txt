############ IO TEST ##############
ADD_EXECUTABLE(load_geo load_geo.cpp ${OpenMEEG_SOURCE_DIR}/src/mesh3.cpp)
TARGET_LINK_LIBRARIES (load_geo ${OPENMEEG_OTHER_LIBRARIES})
ADD_TEST(load_mesh_tri ${EXECUTABLE_OUTPUT_PATH}/load_geo ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri)

ADD_EXECUTABLE(sensors sensors.cpp)
TARGET_LINK_LIBRARIES (sensors ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(compare_matrix compare_matrix.cpp)
TARGET_LINK_LIBRARIES (compare_matrix ${OPENMEEG_OTHER_LIBRARIES})

############ TEST COMMON RESULTS ON HEAD1 (Regression test) ##############

#------ LHS -----#
ADD_TEST(compareLHS ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                    ${OpenMEEG_BINARY_DIR}/tests/lhsmatrix.bin
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_lhsmatrix.bin
                    -sym)
#------ LHS INV -----#
ADD_TEST(compareLHSINV ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                       ${OpenMEEG_BINARY_DIR}/tests/lhsinvmatrix.bin
                       ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_lhsinvmatrix.bin
                       -sym)
#------ RHS -----#
ADD_TEST(compareRHS ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                    ${OpenMEEG_BINARY_DIR}/tests/rhsmatrix.bin
                    ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_rhsmatrix.bin)

############ TEST EEG RESULTS ##############

#------ X2EEG -----#
ADD_TEST(compareX2EEG ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/x2eegmatrix.bin
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_x2eegmatrix.bin)
#------ EEG GAIN -----#
ADD_TEST(compareEEGGAIN ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/eeggainmatrix.bin
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_eeggainmatrix.bin)
#------ EEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareEEGEST ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/estimatedEEG.txt
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedeeg.txt)
#------ EEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareEEGINVHEAT ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                           ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGHEAT.txt
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_eegheat.txt)
#------ EEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareEEGINVMN ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGMN.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_eegmn.txt)
#------ EEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareEEGINVTV ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCEEGTV.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_eegtv.txt)

############ TEST MEG RESULTS ##############

#------ X2MEG -----#
ADD_TEST(compareX2MEG ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/x2megmatrix.bin
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_x2megmatrix.bin)
#------ S2MEG -----#
ADD_TEST(compareS2MEG ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                      ${OpenMEEG_BINARY_DIR}/tests/s2megmatrix.bin
                      ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_s2megmatrix.bin)
#------ MEG GAIN -----#
ADD_TEST(compareMEGGAIN ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/meggainmatrix.bin
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_meggainmatrix.bin)
#------ MEG ESTIMATION (FORWARD) -----#
ADD_TEST(compareMEGEST ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                        ${OpenMEEG_BINARY_DIR}/tests/estimatedMEG.txt
                        ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedmeg.txt)
#------ MEG SOURCE ESTIMATION (INVERSE HEAT) -----#
ADD_TEST(compareMEGINVHEAT ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                           ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGHEAT.txt
                           ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_megheat.txt)
#------ MEG SOURCE ESTIMATION (INVERSE MN) -----#
ADD_TEST(compareMEGINVMN ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGMN.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_megmn.txt)
#------ MEG SOURCE ESTIMATION (INVERSE TV) -----#
ADD_TEST(compareMEGINVTV ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedSRCMEGTV.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/init_estimatedsrc_megtv.txt)

############ TEST EEG RESULTS ON DIPOLES ##############

# Compare EEG result with analytical solution obtained with Matlab
ADD_TEST(compareEEGdipoles ${CMAKE_CURRENT_BINARY_DIR}/compare_matrix
                         ${OpenMEEG_BINARY_DIR}/tests/estimatedEEGpoint.txt
                         ${OpenMEEG_SOURCE_DIR}/tests/initialTest/eeg_analytique_dip.txt
                         -corr -eps 0.005)

