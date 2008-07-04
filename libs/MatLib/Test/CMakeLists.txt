ADD_EXECUTABLE(matlibtest test.cpp chrono.cpp)
TARGET_LINK_LIBRARIES (matlibtest ${OPENMEEG_OTHER_LIBRARIES} OpenMEEGMaths)
ADD_TEST(matlibtest ${CMAKE_CURRENT_BINARY_DIR}/matlibtest)
