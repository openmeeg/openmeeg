include(UnitTest)

OPENMEEG_UNIT_TEST(OpenMEEGMathsTest-vector SOURCES vector.cpp INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} LIBRARIES OpenMEEGMaths)
OPENMEEG_UNIT_TEST(OpenMEEGMathsTest-full SOURCES full.cpp INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} LIBRARIES OpenMEEGMaths)
OPENMEEG_UNIT_TEST(OpenMEEGMathsTest-symm SOURCES symm.cpp INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} LIBRARIES OpenMEEGMaths)
OPENMEEG_UNIT_TEST(OpenMEEGMathsTest-sparse SOURCES sparse.cpp INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} LIBRARIES OpenMEEGMaths)
OPENMEEG_UNIT_TEST(OpenMEEGMathsTest-block SOURCES block.cpp INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} LIBRARIES OpenMEEGMaths)

OPENMEEG_UNIT_TEST(test_mat_files_io
    SOURCES test_mat_files_io.cpp
    LIBRARIES OpenMEEGMaths
    PARAMETERS ${CMAKE_CURRENT_SOURCE_DIR}/data/test_file_v6.mat ${CMAKE_CURRENT_SOURCE_DIR}/data/test_file.mat)

if (NOT matio_FOUND)
    set_tests_properties(test_mat_files_io PROPERTIES WILL_FAIL ON)
endif()
