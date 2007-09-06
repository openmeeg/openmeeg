# A way to load data via their absolute path
# See SRCPATH macro in C++ source
ADD_DEFINITIONS ( -DSRCDIR=${OpenMEEG_SOURCE_DIR} )

ADD_EXECUTABLE(txt2tex txt2tex.cpp)
TARGET_LINK_LIBRARIES (txt2tex ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(tex2txt tex2txt.cpp)
TARGET_LINK_LIBRARIES (tex2txt ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(tex2mask tex2mask.cpp)
TARGET_LINK_LIBRARIES (tex2mask ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(mesh_convert mesh_convert.cpp ../src/mesh3.cpp)
TARGET_LINK_LIBRARIES (mesh_convert ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(mesh_concat mesh_concat.cpp ../src/mesh3.cpp)
TARGET_LINK_LIBRARIES (mesh_concat ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(mattxt2bin mattxt2bin.cpp)
TARGET_LINK_LIBRARIES (mattxt2bin ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(matbin2txt matbin2txt.cpp)
TARGET_LINK_LIBRARIES (matbin2txt ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(symtxt2bin symtxt2bin.cpp)
TARGET_LINK_LIBRARIES (symtxt2bin ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(symbin2txt symbin2txt.cpp)
TARGET_LINK_LIBRARIES (symbin2txt ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(cortex_patch cortex_patch.cpp ../src/mesh3.cpp)
TARGET_LINK_LIBRARIES (cortex_patch ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(register_squids register_squids.cpp)
TARGET_LINK_LIBRARIES (register_squids ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(squids2vtk squids2vtk.cpp)
TARGET_LINK_LIBRARIES (squids2vtk ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(extract_frames extract_frames.cpp)
TARGET_LINK_LIBRARIES (extract_frames ${OPENMEEG_OTHER_LIBRARIES})

ADD_EXECUTABLE(matrix_info matrix_info.cpp)
TARGET_LINK_LIBRARIES (matrix_info ${OPENMEEG_OTHER_LIBRARIES})

IF(ENABLE_PACKAGING)
INSTALL(TARGETS matbin2txt mattxt2bin symtxt2bin symbin2txt mesh_convert mesh_concat cortex_patch register_squids squids2vtk tex2mask tex2txt txt2tex matrix_info DESTINATION bin)
ENDIF(ENABLE_PACKAGING)
