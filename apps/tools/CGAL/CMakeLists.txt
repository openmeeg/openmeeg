if (USE_CGAL)

    find_package(CGAL REQUIRED)

    include_directories(include)

    # Compile the OpenMEEG-CGAL library

    add_library(OpenMEEG-CGAL src/decimate_mesh.cpp src/refine_mesh.cpp src/mesh_from_labeled_image.cpp src/mesh_from_levelset_image.cpp)
    set_target_properties(OpenMEEG-CGAL PROPERTIES VERSION 1.1.0 SOVERSION 1 CLEAN_DIRECT_OUTPUT 1)
    target_link_libraries(OpenMEEG-CGAL OpenMEEG CGAL::CGAL ${VTK_LIBRARIES})
    install(TARGETS OpenMEEG-CGAL DESTINATION lib)

    # Compile each target

    set(OM_CGAL_EXEC)
    foreach(target decimate_mesh refine_mesh mesh_from_levelset_function mesh_from_labeled_image mesh_from_levelset_image)
        add_executable(om_${target} ${target}.cpp)
        list(APPEND OM_CGAL_EXEC om_${target})
        target_link_libraries(om_${target} OpenMEEG OpenMEEG-CGAL CGAL::CGAL ${VTK_LIBRARIES})
    endforeach()
    install(TARGETS ${OM_CGAL_EXEC} DESTINATION bin)

    # tests

    if (BUILD_TESTING)
        OPENMEEG_TEST(cgal_mesh_function om_mesh_from_levelset_function -r 1.0 -fs 10 -fd 1 -o ${OpenMEEG_BINARY_DIR}/tests/sphere_from_levelset.tri)
        OPENMEEG_TEST(cgal_remesh        om_refine_mesh -i ${CMAKE_SOURCE_DIR}/data/Head1/cortex.1.tri -o ${OpenMEEG_BINARY_DIR}/tests/cortex.1.remeshed.tri)
        OPENMEEG_TEST(cgal_decimate      om_decimate_mesh -i ${CMAKE_SOURCE_DIR}/data/Head3/cortex.3.tri -n 100 -o ${OpenMEEG_BINARY_DIR}/tests/cortex.3.decimated.tri)
    endif()
endif()
