IF(NOT WIN32)
    SET(MATIO_EXTERNAL_LIBRARIES m) 
ENDIF()

ADD_EXECUTABLE(test_mat test_mat.c)
TARGET_LINK_LIBRARIES(test_mat matio ${MATIO_EXTERNAL_LIBRARIES})

ADD_EXECUTABLE(test_snprintf test_snprintf.c ${MATIO_SOURCE_DIR}/src/snprintf.c)
TARGET_LINK_LIBRARIES(test_snprintf ${MATIO_EXTERNAL_LIBRARIES})

MACRO(MATIO_TEST NAME PROG_NAME)
    SET(PROG_ARGS ${ARGN})
    IF (WIN32)
        SET (TEST_DIR ${EXECUTABLE_OUTPUT_PATH})
    ELSE ()
        SET (TEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
    ENDIF ()

    SET (EXECUTABLE ${TEST_DIR}/${PROG_NAME}${CMAKE_EXECUTABLE_SUFFIX})
    ADD_TEST(${NAME}-COMPARISON
             ${CMAKE_COMMAND} -D TEST_PROGRAM=${EXECUTABLE}
                              -D TEST_ARGS=${PROG_ARGS}
                              -D TEST_OUTPUT=${TEST_DIR}/${NAME}.out
                              -D TEST_REFERENCE=${MATIO_SOURCE_DIR}/test/results/${NAME}.out
                              -P ${MATIO_SOURCE_DIR}/test/runTest.cmake ) # To compare output to reference file
    ADD_TEST(${NAME} ${EXECUTABLE} ${PROG_ARGS}) # To perform memcheck tests
    #SET_TESTS_PROPERTIES(${NAME} PROPERTIES DEPENDS ${EXECUTABLE})
    SET_TESTS_PROPERTIES(${NAME}-COMPARISON PROPERTIES
                         DEPENDS ${NAME}
                         PASS_REGULAR_EXPRESSION "Success")
ENDMACRO()

MATIO_TEST(TEST_SNPRINTF test_snprintf)

SET(MATIO_WRITE_TESTS write writecompressed writesparse write_compressed_sparse
    write_struct write_compressed_struct writecell write_compressed_cell
    writeinf writenan writeslab)
#    writeinf writenan writenull writeslab)

FOREACH(arg ${MATIO_WRITE_TESTS})
    MATIO_TEST(MATIO-${arg} test_mat ${arg})
ENDFOREACH()

SET(MATIO_IND_TESTS ind2sub sub2ind)
FOREACH(arg ${MATIO_IND_TESTS})
    MATIO_TEST(MATIO-${arg} test_mat ${arg})
ENDFOREACH()

SET(MATIO_FILES write write_compressed writesparse write_compressed_sparse
    write_struct write_compressed_struct writecell write_compressed_cell
    writeinf writenan writeslab)
#    writeinf writenan write_null writeslab)
FOREACH(file ${MATIO_FILES})
    MATIO_TEST(MATIO-copy-${file} test_mat copy test_mat_${file}.mat ${file}.mat)
ENDFOREACH()

SET(MATIO_OTHER_FILES writeinf writenan write_null)
FOREACH(file ${MATIO_FILES})
    MATIO_TEST(MATIO-copy-${file} test_mat copy test_${file}.mat ${file}.mat)
ENDFOREACH()

#   Record the output and compare it with stored values ??

SET(MATIO_WRITE_VARS d f i64 ui64 i32 i16 i8 str)
FOREACH(file write write_compressed)
    FOREACH(var ${MATIO_WRITE_VARS})
        MATIO_TEST(MATIO-readvarinfo-${file}-${var} test_mat readvarinfo test_${file}.mat ${var})
        MATIO_TEST(MATIO-readvar-${file}-${var} test_mat readvar test_${file}.mat ${var})
    ENDFOREACH()
ENDFOREACH()

#FOREACH(file writecell write_compressed_cell)
FOREACH(file writecell)
    MATIO_TEST(MATIO-readvarinfo-${file} test_mat readvarinfo test_mat_${file}.mat cell)
    MATIO_TEST(MATIO-readvar-${file} test_mat readvar test_mat_${file}.mat cell)
ENDFOREACH()

SET(MATIO_WRITESLAB_VARS d f i)
FOREACH(var ${MATIO_WRITESLAB_VARS})
    MATIO_TEST(MATIO-readslab-${var} test_mat readslab test_mat_writeslab.mat ${var})
ENDFOREACH()

FOREACH(file write_struct write_compressed_struct)
    FOREACH(field 1 2 3 4 data)
        MATIO_TEST(MATIO-getstructfield-${file}-${field} test_mat getstructfield test_mat_${file}.mat structure ${var})
    ENDFOREACH()
ENDFOREACH()

# Add more tests for these.

SET(MATIO_READ_TESTS readvar4 readslab4)

# Set tests that are expected to fail (TO BE CORRECTED).

#SET_TESTS_PROPERTIES(MATIO-writenull PROPERTIES WILL_FAIL TRUE)
#SET_TESTS_PROPERTIES(MATIO-writenull-COMPARISON PROPERTIES WILL_FAIL TRUE)
#SET_TESTS_PROPERTIES(MATIO-copy-write_null PROPERTIES WILL_FAIL TRUE)
#SET_TESTS_PROPERTIES(MATIO-readvarinfo-write_compressed_cell PROPERTIES WILL_FAIL TRUE)
#SET_TESTS_PROPERTIES(MATIO-readvar-write_compressed_cell PROPERTIES WILL_FAIL TRUE)

IF (ENABLE_FORTRAN)
    INCLUDE_DIRECTORIES(${MATIO_SOURCE_DIR}/src/fortran ${MATIO_BINARY_DIR}/src/fortran)
    ADD_EXECUTABLE(test_matf test_matf.f90)
    TARGET_LINK_LIBRARIES(test_matf matio m)
    # TESTS
ENDIF()
