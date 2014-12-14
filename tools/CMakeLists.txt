add_executable(om_make_nerve make_nerve.cpp)
target_link_libraries (om_make_nerve OpenMEEG)

add_executable(om_mesh_convert mesh_convert.cpp)
target_link_libraries (om_mesh_convert OpenMEEG)

add_executable(om_geometry_info geometry_info.cpp)
target_link_libraries (om_geometry_info OpenMEEG)

add_executable(om_mesh_info mesh_info.cpp)
target_link_libraries (om_mesh_info OpenMEEG)

add_executable(om_mesh_smooth mesh_smooth.cpp)
target_link_libraries (om_mesh_smooth OpenMEEG)

add_executable(om_mesh_concat mesh_concat.cpp)
target_link_libraries (om_mesh_concat OpenMEEG)

add_executable(om_project_sensors project_sensors.cpp)
target_link_libraries (om_project_sensors OpenMEEG)

add_executable(om_register_squids register_squids.cpp)
target_link_libraries (om_register_squids OpenMEEG)

add_executable(om_squids2vtk squids2vtk.cpp)
target_link_libraries (om_squids2vtk OpenMEEGMaths)

add_executable(om_matrix_info matrix_info.cpp)
target_link_libraries (om_matrix_info OpenMEEGMaths)

add_executable(om_matrix_convert matrix_convert.cpp)
target_link_libraries (om_matrix_convert OpenMEEGMaths)

add_executable(om_check_geom check_geom.cpp)
target_link_libraries (om_check_geom OpenMEEG)

add_executable(om_mesh_to_dip mesh_to_dip.cpp)
target_link_libraries (om_mesh_to_dip OpenMEEG)

include(CgalOption)
if (USE_CGAL)
    add_executable(om_cgal_mesh_3Dlevelset_image cgal_mesh_3Dlevelset_image.cpp)
    target_link_libraries (om_cgal_mesh_3Dlevelset_image OpenMEEG ${CGAL_LIBRARIES})
    add_executable(om_cgal_mesh_3Dlabeled_image cgal_mesh_3Dlabeled_image.cpp)
    target_link_libraries (om_cgal_mesh_3Dlabeled_image OpenMEEG ${CGAL_LIBRARIES})
    add_executable(om_cgal_mesh_function cgal_mesh_function.cpp)
    target_link_libraries (om_cgal_mesh_function OpenMEEG ${CGAL_LIBRARIES})
    add_executable(om_cgal_remesh cgal_remesh.cpp)
    target_link_libraries (om_cgal_remesh OpenMEEG ${CGAL_LIBRARIES})
    install(TARGETS om_cgal_mesh_3Dlevelset_image om_cgal_mesh_3Dlabeled_image om_cgal_mesh_function om_cgal_remesh DESTINATION bin)
endif()

if (USE_VTK)
    add_executable(om_vtp_to_meshes vtp_to_meshes.cpp)
    target_link_libraries (om_vtp_to_meshes OpenMEEG ${VTK_LIBRARIES})

    add_executable(om_meshes_to_vtp meshes_to_vtp.cpp)
    target_link_libraries (om_meshes_to_vtp OpenMEEG ${VTK_LIBRARIES})

    add_executable(om_add_dataset_to_vtk add_dataset_to_vtk.cpp)
    target_link_libraries (om_add_dataset_to_vtk OpenMEEG ${VTK_LIBRARIES})
    install(TARGETS om_vtp_to_meshes om_meshes_to_vtp om_add_dataset_to_vtk DESTINATION bin)
endif()

if (ENABLE_PYTHON)
    subdirs(python)
endif()

install(TARGETS om_make_nerve om_mesh_convert om_mesh_concat om_project_sensors om_mesh_info om_mesh_smooth om_register_squids om_geometry_info om_squids2vtk om_matrix_info om_matrix_convert om_check_geom om_mesh_to_dip DESTINATION bin)

if (USE_VTK)
    OPENMEEG_TEST(Tool-om_meshes_to_vtp ${CMAKE_CURRENT_BINARY_DIR}/om_meshes_to_vtp -i1 ${CMAKE_SOURCE_DIR}/data/Head1/cortex.1.tri -i2 ${CMAKE_SOURCE_DIR}/data/Head1/skull.1.tri -i3 ${CMAKE_SOURCE_DIR}/data/Head1/scalp.1.tri -n1 "cortex" -n2 "skull" -n3 "scalp" -o ${OpenMEEG_BINARY_DIR}/tests/Head1.vtp)
    OPENMEEG_TEST(Tool-om_vtp_to_meshes ${CMAKE_CURRENT_BINARY_DIR}/om_vtp_to_meshes -i ${OpenMEEG_BINARY_DIR}/tests/Head1.vtp -o ${OpenMEEG_BINARY_DIR}/tests/meshes.tri)
endif()
