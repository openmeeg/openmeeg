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

if (USE_CGAL)
    set(OM_CGAL_TARGETS cgal_mesh_function)
    set(OM_CGAL_SOURCES cgal_lib_function.cpp)
    include_directories(.)

    if (CGAL_ImageIO_FOUND) # these are tools handling images (.inr format only for the moment!)
        list(APPEND OM_CGAL_TARGETS "cgal_mesh_3Dlabeled_image" "cgal_mesh_3Dlevelset_image")
        list(APPEND OM_CGAL_EXEC om_${target})
        list(APPEND OM_CGAL_SOURCES cgal_lib_image.cpp)
    endif()

    if (CGAL_VERSION VERSION_GREATER 4.6)
        # edge collapse feature is stabilized in version ? not 4.2
        list(APPEND OM_CGAL_TARGETS "cgal_decimate" "cgal_remesh")
        list(APPEND OM_CGAL_SOURCES cgal_lib_decimate.cpp)
    endif()

    # compile the OpenMEEGCgal library
    add_library(OpenMEEGCgal ${OM_CGAL_SOURCES})

    set_target_properties(OpenMEEGCgal PROPERTIES
                          VERSION 1.1.0
                          SOVERSION 1
                          CLEAN_DIRECT_OUTPUT 1)

    target_link_libraries(OpenMEEGCgal OpenMEEG ${CGAL_LIBRARIES} ${VTK_LIBRARIES})
    target_include_directories(OpenMEEGCgal PUBLIC ${CGAL_INCLUDE_DIRS})
    target_compile_options(OpenMEEGCgal PUBLIC ${CGAL_CXX_FLAGS})
    install(TARGETS OpenMEEGCgal DESTINATION lib)

    # compile each target
    set(OM_CGAL_EXEC)
    foreach(target ${OM_CGAL_TARGETS})
        add_executable(om_${target} ${target}.cpp)
        list(APPEND OM_CGAL_EXEC om_${target})
        target_link_libraries(om_${target} OpenMEEG OpenMEEGCgal ${CGAL_LIBRARIES} ${VTK_LIBRARIES})
        target_include_directories(om_${target} PUBLIC ${CGAL_INCLUDE_DIRS})
        target_compile_options(om_${target} PUBLIC ${CGAL_CXX_FLAGS})
    endforeach()
    install(TARGETS ${OM_CGAL_EXEC} DESTINATION bin)

    # tests
    if (BUILD_TESTING)
        OPENMEEG_TEST(cgal_mesh_function om_cgal_mesh_function -r 1.0 -fs 10 -fd 1 -o ${OpenMEEG_BINARY_DIR}/tests/remeshed.tri)
        if (CGAL_VERSION VERSION_GREATER 4.6)
            OPENMEEG_TEST(cgal_remesh om_cgal_remesh -i ${CMAKE_SOURCE_DIR}/data/Head1/cortex.1.tri -o ${OpenMEEG_BINARY_DIR}/tests/remeshed.tri)
            OPENMEEG_TEST(cgal_decimate om_cgal_decimate -i ${CMAKE_SOURCE_DIR}/data/Head3/cortex.3.tri -n 100 -o ${OpenMEEG_BINARY_DIR}/tests/decimated.tri)
        endif()
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

# tests
if (BUILD_TESTING)
    OPENMEEG_TEST(wrong_geom_info om_geometry_info -g ${CMAKE_SOURCE_DIR}/data/HeadNNc1/HeadNNc1.geom -c ${CMAKE_SOURCE_DIR}/data/Head1/Head1.cond)
    set_tests_properties(wrong_geom_info PROPERTIES WILL_FAIL TRUE)
endif()

install(TARGETS om_make_nerve om_mesh_convert om_mesh_concat om_project_sensors om_mesh_info om_mesh_smooth om_register_squids om_geometry_info om_squids2vtk om_matrix_info om_matrix_convert om_check_geom om_mesh_to_dip DESTINATION bin)
