ADD_EXECUTABLE(om_make_nerve make_nerve.cpp)
TARGET_LINK_LIBRARIES (om_make_nerve OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_mesh_convert mesh_convert.cpp)
TARGET_LINK_LIBRARIES (om_mesh_convert OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_mesh_info mesh_info.cpp)
TARGET_LINK_LIBRARIES (om_mesh_info OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_mesh_smooth mesh_smooth.cpp)
TARGET_LINK_LIBRARIES (om_mesh_smooth OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_mesh_concat mesh_concat.cpp)
TARGET_LINK_LIBRARIES (om_mesh_concat OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_register_squids register_squids.cpp)
TARGET_LINK_LIBRARIES (om_register_squids OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_squids2vtk squids2vtk.cpp)
TARGET_LINK_LIBRARIES (om_squids2vtk OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_matrix_info matrix_info.cpp)
TARGET_LINK_LIBRARIES (om_matrix_info OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_matrix_convert matrix_convert.cpp)
TARGET_LINK_LIBRARIES (om_matrix_convert OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_vector_convert vector_convert.cpp)
TARGET_LINK_LIBRARIES (om_vector_convert OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_symmatrix_convert symmatrix_convert.cpp)
TARGET_LINK_LIBRARIES (om_symmatrix_convert OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_sparse_matrix_convert sparse_matrix_convert.cpp)
TARGET_LINK_LIBRARIES (om_sparse_matrix_convert OpenMEEGMaths "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_check_geom check_geom.cpp)
TARGET_LINK_LIBRARIES (om_check_geom OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

ADD_EXECUTABLE(om_mesh_to_dip mesh_to_dip.cpp)
TARGET_LINK_LIBRARIES (om_mesh_to_dip OpenMEEG "${OPENMEEG_OTHER_LIBRARIES}")

INSTALL(TARGETS om_make_nerve om_mesh_convert om_mesh_concat om_mesh_info om_mesh_smooth om_register_squids om_squids2vtk om_matrix_info om_matrix_convert om_symmatrix_convert om_sparse_matrix_convert om_vector_convert om_check_geom om_mesh_to_dip DESTINATION bin)
