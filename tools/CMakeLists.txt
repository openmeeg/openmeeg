# A way to load data via their absolute path
# See SRCPATH macro in C++ source

IF ( NOT XCODE )
    ADD_DEFINITIONS ( -D"SRCDIR=${OpenMEEG_SOURCE_DIR}" )
ENDIF ( NOT XCODE )

ADD_EXECUTABLE(om_txt2tex txt2tex.cpp)
TARGET_LINK_LIBRARIES (om_txt2tex ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_make_nerve make_nerve.cpp)
TARGET_LINK_LIBRARIES (om_make_nerve OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_tex2txt tex2txt.cpp)
TARGET_LINK_LIBRARIES (om_tex2txt ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_tex2mask tex2mask.cpp)
TARGET_LINK_LIBRARIES (om_tex2mask ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_mesh_convert mesh_convert.cpp)
TARGET_LINK_LIBRARIES (om_mesh_convert OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_mesh_info mesh_info.cpp)
TARGET_LINK_LIBRARIES (om_mesh_info OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_mesh_smooth mesh_smooth.cpp)
TARGET_LINK_LIBRARIES (om_mesh_smooth OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_mesh_concat mesh_concat.cpp)
TARGET_LINK_LIBRARIES (om_mesh_concat OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_mattxt2bin mattxt2bin.cpp)
TARGET_LINK_LIBRARIES (om_mattxt2bin ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_matbin2txt matbin2txt.cpp)
TARGET_LINK_LIBRARIES (om_matbin2txt ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_symtxt2bin symtxt2bin.cpp)
TARGET_LINK_LIBRARIES (om_symtxt2bin ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_symbin2txt symbin2txt.cpp)
TARGET_LINK_LIBRARIES (om_symbin2txt ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_cortex_patch cortex_patch.cpp)
TARGET_LINK_LIBRARIES (om_cortex_patch OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_register_squids register_squids.cpp)
TARGET_LINK_LIBRARIES (om_register_squids OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_squids2vtk squids2vtk.cpp)
TARGET_LINK_LIBRARIES (om_squids2vtk ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_extract_frames extract_frames.cpp)
TARGET_LINK_LIBRARIES (om_extract_frames ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_matrix_info matrix_info.cpp)
TARGET_LINK_LIBRARIES (om_matrix_info ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_matrix_convert matrix_convert.cpp)
TARGET_LINK_LIBRARIES (om_matrix_convert ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_check_geom check_geom.cpp)
TARGET_LINK_LIBRARIES (om_check_geom OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(om_mesh_to_dip mesh_to_dip.cpp)
TARGET_LINK_LIBRARIES (om_mesh_to_dip OpenMEEG ${OPENMEEG_OTHER_LIBRARIES})

IF(ENABLE_PACKAGING)
INSTALL(TARGETS om_make_nerve om_matbin2txt om_mattxt2bin om_symtxt2bin om_symbin2txt om_mesh_convert om_mesh_concat om_mesh_info om_mesh_smooth om_cortex_patch om_register_squids om_squids2vtk om_tex2mask om_tex2txt om_txt2tex om_matrix_info om_matrix_convert om_check_geom om_mesh_to_dip DESTINATION bin)
ENDIF(ENABLE_PACKAGING)
