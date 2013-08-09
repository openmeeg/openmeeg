#ADD_EXECUTABLE(om_make_nerve make_nerve.cpp)
#TARGET_LINK_LIBRARIES (om_make_nerve OpenMEEG)

ADD_EXECUTABLE(om_mesh_convert mesh_convert.cpp)
TARGET_LINK_LIBRARIES (om_mesh_convert OpenMEEG)

ADD_EXECUTABLE(om_geometry_info geometry_info.cpp)
TARGET_LINK_LIBRARIES (om_geometry_info OpenMEEG)

ADD_EXECUTABLE(om_mesh_info mesh_info.cpp)
TARGET_LINK_LIBRARIES (om_mesh_info OpenMEEG)

ADD_EXECUTABLE(om_mesh_smooth mesh_smooth.cpp)
TARGET_LINK_LIBRARIES (om_mesh_smooth OpenMEEG)

#ADD_EXECUTABLE(om_mesh_concat mesh_concat.cpp)
#TARGET_LINK_LIBRARIES (om_mesh_concat OpenMEEG)

ADD_EXECUTABLE(om_register_squids register_squids.cpp)
TARGET_LINK_LIBRARIES (om_register_squids OpenMEEG)

ADD_EXECUTABLE(om_squids2vtk squids2vtk.cpp)
TARGET_LINK_LIBRARIES (om_squids2vtk OpenMEEGMaths)

ADD_EXECUTABLE(om_matrix_info matrix_info.cpp)
TARGET_LINK_LIBRARIES (om_matrix_info OpenMEEGMaths)

ADD_EXECUTABLE(om_matrix_convert matrix_convert.cpp)
TARGET_LINK_LIBRARIES (om_matrix_convert OpenMEEGMaths)

ADD_EXECUTABLE(om_check_geom check_geom.cpp)
TARGET_LINK_LIBRARIES (om_check_geom OpenMEEG)

ADD_EXECUTABLE(om_vtp_to_mesh vtp_to_mesh.cpp)
TARGET_LINK_LIBRARIES (om_vtp_to_mesh OpenMEEG)

ADD_EXECUTABLE(om_mesh_to_vtp mesh_to_vtp.cpp)
TARGET_LINK_LIBRARIES (om_mesh_to_vtp OpenMEEG)

ADD_EXECUTABLE(om_mesh_to_dip mesh_to_dip.cpp)
TARGET_LINK_LIBRARIES (om_mesh_to_dip OpenMEEG)

INCLUDE(CgalOption)
IF (USE_CGAL)
    ADD_EXECUTABLE(om_cgal_mesh_create cgal_mesh_create.cpp)
    TARGET_LINK_LIBRARIES (om_cgal_mesh_create OpenMEEG ${CGAL_LIBRARIES})
    INSTALL(TARGETS om_cgal_mesh_create DESTINATION bin)
ENDIF()

IF (USE_VTK)
    ADD_EXECUTABLE(om_add_vtk_dataset AddDataSetToVtk.cpp)
    TARGET_LINK_LIBRARIES (om_add_vtk_dataset ${VTK_LIBRARIES})
    INSTALL(TARGETS om_add_vtk_dataset DESTINATION bin)
ENDIF()

INSTALL(TARGETS 
    #om_make_nerve
    om_mesh_convert 
    # om_mesh_concat
    om_mesh_info om_mesh_smooth om_register_squids
    om_geometry_info
    om_mesh_to_vtp om_vtp_to_mesh
    om_squids2vtk om_matrix_info om_matrix_convert om_check_geom om_mesh_to_dip DESTINATION bin
    )

OPENMEEG_TEST(Tool-om_mesh_to_vtp ${CMAKE_CURRENT_BINARY_DIR}/om_mesh_to_vtp -i1 ${CMAKE_SOURCE_DIR}/data/Models/Head1/cortex.1.tri -i2 ${CMAKE_SOURCE_DIR}/data/Models/Head1/skull.1.tri -i3 ${CMAKE_SOURCE_DIR}/data/Models/Head1/scalp.1.tri -n1 "cortex" -n2 "skull" -n3 "scalp" -o ${CMAKE_CURRENT_BINARY_DIR}/Head1.vtp)

OPENMEEG_TEST(Tool-om_vtp_to_mesh ${CMAKE_CURRENT_BINARY_DIR}/om_vtp_to_mesh -i ${CMAKE_CURRENT_BINARY_DIR}/Head1.vtp -o ${CMAKE_CURRENT_BINARY_DIR}/meshes.tri)
