if (USE_VTK)
    add_executable(om_vtp_to_meshes vtp_to_meshes.cpp)
    target_link_libraries(om_vtp_to_meshes OpenMEEG ${VTK_LIBRARIES})

    add_executable(om_meshes_to_vtp meshes_to_vtp.cpp)
    target_link_libraries(om_meshes_to_vtp OpenMEEG ${VTK_LIBRARIES})

    add_executable(om_add_dataset_to_vtk add_dataset_to_vtk.cpp)
    target_link_libraries(om_add_dataset_to_vtk OpenMEEG ${VTK_LIBRARIES})

    install(TARGETS om_vtp_to_meshes om_meshes_to_vtp om_add_dataset_to_vtk DESTINATION bin)

    if (BUILD_TESTING)
        OPENMEEG_TEST(meshes_to_vtp om_meshes_to_vtp -g ${CMAKE_SOURCE_DIR}/data/Head1/Head1.geom -o ${OpenMEEG_BINARY_DIR}/tests/Head1.vtp)
        OPENMEEG_TEST(vtp_to_meshes om_vtp_to_meshes -i ${OpenMEEG_BINARY_DIR}/tests/Head1.vtp -o ${OpenMEEG_BINARY_DIR}/tests/meshes.tri DEPENDS meshes_to_vtp)
    endif()
endif()
