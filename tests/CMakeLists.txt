############ IO TEST ##############
ADD_EXECUTABLE(load_geo load_geo.cpp ${OpenMEEG_SOURCE_DIR}/src/mesh3.cpp)
TARGET_LINK_LIBRARIES (load_geo ${OPENMEEG_OTHER_LIBRARIES})
ADD_TEST(load_tri ${EXECUTABLE_OUTPUT_PATH}/load_geo ${OpenMEEG_SOURCE_DIR}/HeadModels/Head1/Head1.tri)
