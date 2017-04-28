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
target_link_libraries(om_squids2vtk OpenMEEGMaths)

add_executable(om_matrix_info matrix_info.cpp)
target_link_libraries(om_matrix_info OpenMEEGMaths)

add_executable(om_matrix_convert matrix_convert.cpp)
target_link_libraries(om_matrix_convert OpenMEEGMaths)

add_executable(om_check_geom check_geom.cpp)
target_link_libraries(om_check_geom OpenMEEG ${VTK_LIBRARIES})

add_executable(om_mesh_to_dip mesh_to_dip.cpp)
target_link_libraries(om_mesh_to_dip OpenMEEG ${VTK_LIBRARIES})

if (USE_CGAL)
    set(OM_CGAL_EXEC)
    foreach(target cgal_mesh_function cgal_remesh)
        add_executable(om_${target} ${target}.cpp)
        list(APPEND OM_CGAL_EXEC om_${target})
        target_link_libraries(om_${target} OpenMEEG ${CGAL_LIBRARIES} ${VTK_LIBRARIES})
        target_include_directories(om_${target} PUBLIC ${CGAL_INCLUDE_DIRS})
        target_compile_options(om_${target} PUBLIC ${CGAL_CXX_FLAGS})
    endforeach()
    if (CGAL_ImageIO_FOUND) # these are tools hanling images (.inr format only for the moment!)
        foreach(target cgal_mesh_3Dlevelset_image cgal_mesh_3Dlabeled_image)
            add_executable(om_${target} ${target}.cpp)
            list(APPEND OM_CGAL_EXEC om_${target})
            target_link_libraries(om_${target} OpenMEEG ${CGAL_LIBRARIES} ${VTK_LIBRARIES})
            target_include_directories(om_${target} PUBLIC ${CGAL_INCLUDE_DIRS})
            target_compile_options(om_${target} PUBLIC ${CGAL_CXX_FLAGS})
        endforeach()
    endif()
    install(TARGETS ${OM_CGAL_EXEC} DESTINATION bin)
    if (BUILD_TESTING)
        OPENMEEG_TEST(cgal_remesh om_cgal_remesh -i ${CMAKE_SOURCE_DIR}/data/Head1/cortex.1.tri -o ${OpenMEEG_BINARY_DIR}/tests/remeshed.tri)
    endif()
endif()

if (USE_VTK)
    add_executable(om_vtp_to_meshes vtp_to_meshes.cpp)
    target_link_libraries(om_vtp_to_meshes OpenMEEG ${VTK_LIBRARIES})

    add_executable(om_meshes_to_vtp meshes_to_vtp.cpp)
    target_link_libraries(om_meshes_to_vtp OpenMEEG ${VTK_LIBRARIES})

    add_executable(om_add_dataset_to_vtk add_dataset_to_vtk.cpp)
    target_link_libraries(om_add_dataset_to_vtk OpenMEEG ${VTK_LIBRARIES})

    install(TARGETS om_vtp_to_meshes om_meshes_to_vtp om_add_dataset_to_vtk 
            DESTINATION bin)

    if (BUILD_TESTING)
        OPENMEEG_TEST(meshes_to_vtp om_meshes_to_vtp -g ${CMAKE_SOURCE_DIR}/data/Head1/Head1.geom -o ${OpenMEEG_BINARY_DIR}/tests/Head1.vtp)
        OPENMEEG_TEST(vtp_to_meshes om_vtp_to_meshes -i ${OpenMEEG_BINARY_DIR}/tests/Head1.vtp -o ${OpenMEEG_BINARY_DIR}/tests/meshes.tri DEPENDS meshes_to_vtp)
    endif()
endif()

# test mesh IOs
if (BUILD_TESTING)
    # extensions to test
    set(extensions off bnd mesh)
    if (USE_VTK)
        list(APPEND extensions vtk)
    endif()
    if (USE_GIFTI)
        list(APPEND extensions gii)
    endif()
    # from tri -> to ext
    foreach(extension ${extensions})
        OPENMEEG_TEST(mesh_convert_tri_to_${extension} om_mesh_convert -i ${CMAKE_SOURCE_DIR}/data/Head1/scalp.1.tri -o ${OpenMEEG_BINARY_DIR}/tests/scalp.1.${extension})
    endforeach()
    # from ext -> to tri
    foreach(extension ${extensions})
        OPENMEEG_TEST(mesh_convert_${extension}_to_tri om_mesh_convert -i ${OpenMEEG_BINARY_DIR}/tests/scalp.1.${extension} -o ${OpenMEEG_BINARY_DIR}/tests/scalp.1.tri)
    endforeach()
endif()

install(TARGETS om_make_nerve om_mesh_convert om_mesh_concat om_project_sensors om_mesh_info om_mesh_smooth om_register_squids om_geometry_info om_squids2vtk om_matrix_info om_matrix_convert om_check_geom om_mesh_to_dip DESTINATION bin)
