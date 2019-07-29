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

# tests
if (BUILD_TESTING)
    OPENMEEG_TEST(check_test_load_geo
        test_load_geo ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.cond)
    OPENMEEG_TEST(check_test_mesh_ios
        test_mesh_ios ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.tri)
endif()


######

