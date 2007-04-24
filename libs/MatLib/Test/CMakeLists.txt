
ADD_EXECUTABLE(matlibtest test.cpp chrono.cpp)
TARGET_LINK_LIBRARIES (matlibtest ${OPENMEEG_OTHER_LIBRARIES})
ADD_TEST(matlibtest ${EXECUTABLE_OUTPUT_PATH}/matlibtest)

ADD_EXECUTABLE(comparematrix compareMatrix.cpp)
TARGET_LINK_LIBRARIES (comparematrix ${OPENMEEG_OTHER_LIBRARIES})

