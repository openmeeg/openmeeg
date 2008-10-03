ADD_EXECUTABLE(matlibtest test.cpp chrono.cpp)
TARGET_LINK_LIBRARIES (matlibtest ${OPENMEEG_OTHER_LIBRARIES} OpenMEEGMaths)

IF (WIN32)
    ADD_TEST(matlibtest ${EXECUTABLE_OUTPUT_PATH}/matlibtest)
ELSE (WIN32)
    ADD_TEST(matlibtest ${CMAKE_CURRENT_BINARY_DIR}/matlibtest)
ENDIF (WIN32)
