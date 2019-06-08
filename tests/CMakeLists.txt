# OPENMEEG_UNIT_TEST(load_geo
#     SOURCES load_geo.cpp
#     LIBRARIES OpenMEEG OpenMEEGMaths ${VTK_LIBRARIES}
#     PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.geom ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.cond)
# 
# OPENMEEG_UNIT_TEST(test_mesh_ios
#     SOURCES test_mesh_ios.cpp
#     LIBRARIES OpenMEEG OpenMEEGMaths ${VTK_LIBRARIES}
#     PARAMETERS ${OpenMEEG_SOURCE_DIR}/data/Head1/Head1.tri)
# 
# #   Test sensor is not used  Why ?
# 
# NEW_EXECUTABLE(test_sensors test_sensors.cpp LIBRARIES OpenMEEG ${VTK_LIBRARIES})
# NEW_EXECUTABLE(compare_matrix compare_matrix.cpp LIBRARIES OpenMEEGMaths ${LAPACK_LIBRARIES})
# NEW_EXECUTABLE(om_validationEIT validationEIT.cpp ${OPEMEEG_HEADERS} LIBRARIES OpenMEEG OpenMEEGMaths ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(load_geo load_geo.cpp)
target_link_libraries(load_geo OpenMEEG::OpenMEEG)

add_executable(test_mesh_ios test_mesh_ios.cpp)
target_link_libraries(test_mesh_ios OpenMEEG::OpenMEEG)

add_executable(test_sensors test_sensors.cpp)
target_link_libraries(test_sensors OpenMEEG::OpenMEEG)

add_executable(compare_matrix compare_matrix.cpp)
target_link_libraries(compare_matrix OpenMEEG::OpenMEEG OpenMEEG::OpenMEEGMaths)

add_executable(om_validationEIT validationEIT.cpp)
target_link_libraries(om_validationEIT OpenMEEG::OpenMEEG)
