add_executable(test_load_geo test_load_geo.cpp)
target_link_libraries(test_load_geo OpenMEEG::OpenMEEG)

add_executable(test_mesh_ios test_mesh_ios.cpp)
target_link_libraries(test_mesh_ios OpenMEEG::OpenMEEG)

add_executable(test_sensors test_sensors.cpp)
target_link_libraries(test_sensors OpenMEEG::OpenMEEG)

add_executable(test_validationEIT test_validationEIT.cpp)
target_link_libraries(test_validationEIT OpenMEEG::OpenMEEG)

add_executable(test_compare_matrix test_compare_matrix.cpp)
target_link_libraries(test_compare_matrix OpenMEEG::OpenMEEG OpenMEEG::OpenMEEGMaths)

OPENMEEG_TEST(check_test_load_geo_legacy
    test_load_geo ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1_legacy.geom ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.cond)
OPENMEEG_TEST(check_test_load_geo
    test_load_geo ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.cond)
OPENMEEG_TEST(check_test_mesh_ios
    test_mesh_ios ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.tri)

include(TestHead.cmake)

# Those models should give same results !
#TESTHEAD(0)
#TESTHEAD(NNa0)
TESTHEAD(1)
TESTHEAD(NNa1)
TESTHEAD(NNb1)
TESTHEAD(NNc1)

TESTHEAD(2)
if (USE_VTK)
    TESTHEAD(NNa2)
    TESTHEAD(NNb2)
    TESTHEAD(NNc2)
endif()

if (TEST_HEAD3)
    TESTHEAD(3)
    if (USE_VTK)
        TESTHEAD(NNa3)
        TESTHEAD(NNb3)
        TESTHEAD(NNc3)
    endif()
endif()

# Multiple nonconductive tests

TESTHEAD(MN1)
TESTHEAD(MN2)
TESTHEAD(MN3)

include(ComparisonTest_caller.cmake)
