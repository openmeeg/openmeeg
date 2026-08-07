add_executable(om_make_nerve make_nerve.cpp)
target_link_libraries(om_make_nerve OpenMEEG ${VTK_LIBRARIES})

add_executable(om_mesh_convert mesh_convert.cpp)
target_link_libraries(om_mesh_convert OpenMEEG ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(om_geometry_info geometry_info.cpp)
target_link_libraries(om_geometry_info OpenMEEG ${VTK_LIBRARIES})

add_executable(om_mesh_info mesh_info.cpp)
target_link_libraries(om_mesh_info OpenMEEG ${VTK_LIBRARIES})

add_executable(om_mesh_smooth mesh_smooth.cpp)
target_link_libraries(om_mesh_smooth OpenMEEG ${VTK_LIBRARIES})

add_executable(om_mesh_concat mesh_concat.cpp)
target_link_libraries(om_mesh_concat OpenMEEG ${VTK_LIBRARIES})

add_executable(om_project_sensors project_sensors.cpp)
target_link_libraries(om_project_sensors OpenMEEG ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(om_register_squids register_squids.cpp)
target_link_libraries(om_register_squids OpenMEEG ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(om_squids2vtk squids2vtk.cpp)
target_link_libraries(om_squids2vtk OpenMEEGMaths OpenMEEG::OpenMEEG)

add_executable(om_matrix_info matrix_info.cpp)
target_link_libraries(om_matrix_info OpenMEEGMaths OpenMEEG::OpenMEEG)

add_executable(om_matrix_convert matrix_convert.cpp)
target_link_libraries(om_matrix_convert OpenMEEGMaths OpenMEEG::OpenMEEG)

add_executable(om_check_geom check_geom.cpp)
target_link_libraries(om_check_geom OpenMEEG ${VTK_LIBRARIES})

# tests
if (BUILD_TESTING)
    OPENMEEG_TEST(check_geom_mesh
        om_check_geom -g ${CMAKE_SOURCE_DIR}/data/Head1/Head1.geom -m ${CMAKE_SOURCE_DIR}/data/Head1/Head1.tri)
    OPENMEEG_TEST(check_geom_dip
        om_check_geom -g ${CMAKE_SOURCE_DIR}/data/Head1/Head1.geom -d ${CMAKE_SOURCE_DIR}/data/Head1/Head1.dip)
    OPENMEEG_TEST(check_geom_dip_outside
        om_check_geom -g ${CMAKE_SOURCE_DIR}/data/Head1/Head1.geom -d ${CMAKE_SOURCE_DIR}/data/Head1/Head1_outside.dip)
    set_tests_properties(check_geom_dip_outside PROPERTIES WILL_FAIL TRUE)
endif()

add_executable(om_mesh_to_dip mesh_to_dip.cpp)
target_link_libraries(om_mesh_to_dip OpenMEEG ${VTK_LIBRARIES})

add_subdirectory(CGAL)
add_subdirectory(vtk)

# tests
if (BUILD_TESTING)
    OPENMEEG_TEST(wrong_geom_info om_geometry_info -g ${CMAKE_SOURCE_DIR}/data/HeadNNc1/HeadNNc1.geom -c ${CMAKE_SOURCE_DIR}/data/Head1/Head1.cond)
    set_tests_properties(wrong_geom_info PROPERTIES WILL_FAIL TRUE)
endif()

install(TARGETS om_make_nerve om_mesh_convert om_mesh_concat om_project_sensors om_mesh_info om_mesh_smooth om_register_squids om_geometry_info om_squids2vtk om_matrix_info om_matrix_convert om_check_geom om_mesh_to_dip DESTINATION bin)
