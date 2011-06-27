IF(NOT WIN32)
    SET(MATIO_EXTERNAL_LIBRARIES m) 
ENDIF()

ADD_EXECUTABLE(test_mat test_mat.c)
TARGET_LINK_LIBRARIES(test_mat matio ${MATIO_EXTERNAL_LIBRARIES} ${GETOPT_LIB})

ADD_EXECUTABLE(test_snprintf test_snprintf.c ${MATIO_SOURCE_DIR}/src/snprintf.c)
TARGET_LINK_LIBRARIES(test_snprintf ${MATIO_EXTERNAL_LIBRARIES} ${GETOPT_LIB})

MACRO(PARSE_TEST_ARGUMENTS LIST_VARS DEFAULT_VAR)
    foreach(var ${LIST_VARS})
        unset(${var})
    endforeach ()

    set(CURRENT_VAR ${DEFAULT_VAR})
    foreach (arg ${ARGN})
        set(skip_this_arg FALSE)
        foreach(var ${LIST_VARS})
            if (${arg} STREQUAL ${var})
                set(CURRENT_VAR ${var})
                set(skip_this_arg TRUE)
                BREAK()
            endif()
        endforeach ()
        if (NOT skip_this_arg)
            set(${CURRENT_VAR} ${${CURRENT_VAR}} ${arg})
        endif()
    endforeach ()
ENDMACRO()

MACRO(MATIO_TEST NAME FILENAME OUTPUT PROG_NAME)
    PARSE_TEST_ARGUMENTS("COMMAND;DEPENDS" "DEFAULT" ${ARGN})
    SET(PROG_ARGS "${ARGN}")
    SEPARATE_ARGUMENTS(ARGS UNIX_COMMAND "${PROG_ARGS}")
    IF (WIN32)
        IF (NOT CMAKE_BUILD_TYPE)
            SET(CMAKE_BUILD_TYPE Debug)
        ENDIF()
        SET (TEST_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_BUILD_TYPE})
    ELSE ()
        SET (TEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
    ENDIF ()

    SET (EXECUTABLE ${TEST_DIR}/${PROG_NAME}${CMAKE_EXECUTABLE_SUFFIX})
    ADD_TEST(${NAME} ${EXECUTABLE} ${PROG_ARGS} -o ${FILENAME}) # To perform memcheck tests
    ADD_TEST(${NAME}-COMPARISON
             ${CMAKE_COMMAND} -D TEST_DIR:STRING=${TEST_DIR}
                              -D TEST_NAME:STRING=${NAME}
                              -D TEST_REFERENCE_DIR:STRING=${MATIO_SOURCE_DIR}/test/results
                              -D TEST_PROGRAM:STRING=${EXECUTABLE}
                              -D TEST_ARGS:STRING=${ARGS}
                              -D TEST_RESULT:STRING=${FILENAME}
                              -D OUTPUT:STRING=${OUTPUT}
                              -P ${MATIO_SOURCE_DIR}/test/runTest.cmake) # To compare output to reference file
    #SET_TESTS_PROPERTIES(${NAME} PROPERTIES DEPENDS ${EXECUTABLE})
    SET_TESTS_PROPERTIES(${NAME}-COMPARISON PROPERTIES
                         DEPENDS ${NAME}
                         PASS_REGULAR_EXPRESSION "Success")
ENDMACRO()

#MATIO_TEST(TEST_SNPRINTF test_snprintf)

SET(MATIO_WRITE_TESTS_NUMERIC 
     write_2d_numeric write_complex_2d_numeric write_empty_2d_numeric write_char write_struct_2d_numeric
     write_struct_complex_2d_numeric write_cell_2d_numeric write_cell_complex_2d_numeric write_char)

SET(MATIO_WRITE_TESTS_SPARSE write_sparse write_complex_sparse)

SET(MATIO_WRITE_TESTS ${MATIO_WRITE_TESTS_NUMERIC} ${MATIO_WRITE_TESTS_SPARSE}
    write_char write_empty_struct 
    write_empty_cell
    writeinf writenull writenan writeslab)

SET(MATIO_FILES)
FOREACH(arg ${MATIO_WRITE_TESTS})
    SET(VERSIONS 4 5)
    IF (NOT ${arg} STREQUAL "writeslab")
        IF (MAT73)
            SET(VERSIONS ${VERSIONS} 7.3)
        ENDIF()
    ENDIF()
    FOREACH(class double single int64 uint64 int32 uint32 int16 uint16 int8 uint8)
        FOREACH (version ${VERSIONS})
            SET(testname ${arg}-${class}-${version})
            SET(filename test_${testname}.mat)
            SET(MATIO_FILES ${MATIO_FILES} ${filename})
            MATIO_TEST(MATIO-${testname} ${filename} Null.out test_mat -c ${class} -v ${version} ${arg} -o ${filename})
            IF (${version} STREQUAL "5")
                SET(filename test_${testname}-compressed.mat)
                MATIO_TEST(MATIO-${testname}-compressed ${filename} Null.out test_mat -c ${class} -v ${version} -z ${arg} -o ${filename})
            ENDIF() 
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

SET(MATIO_OTHER_FILES writeinf writenan writenull writeslab)
FOREACH(arg ${MATIO_OTHER_FILES})
    SET(testname ${arg})
    SET(filename test_${testname}.mat)
    IF (${arg} STREQUAL "writenull")
        SET(filename test_write_null.mat)
    ENDIF()
    SET(MATIO_FILES ${MATIO_FILES} ${filename})
    MATIO_TEST(MATIO-${testname} ${filename} Null.out test_mat ${arg})
ENDFOREACH()

SET(MATIO_IND_TESTS ind2sub sub2ind)
FOREACH(arg ${MATIO_IND_TESTS})
    SET(TEST_REFERENCE MATIO-${arg}.out)
    IF (WIN32)
        SET(TEST_REFERENCE MATIO-${arg}-win.out)
    ENDIF()
    MATIO_TEST(MATIO-${arg} NONE ${TEST_REFERENCE} test_mat ${arg})
ENDFOREACH()

#SET(DATASETS d f i64 ui64 i32 i16 i8 str)
#FOREACH(file ${MATIO_FILES})
#    MATIO_TEST(MATIO-copy-${file} copy_${file} Copy.out test_mat copy ${file} -o copy_${file})
#    FOREACH (var ${DATASETS})
#        ADD_TEST(MATIO-delete-${file} test_mat delete ${file} ${var})
#    ENDFOREACH()
#ENDFOREACH()

#SET(MATIO_WRITESLAB_VARS d f i)
#FOREACH(var ${MATIO_WRITESLAB_VARS})
#    MATIO_TEST(MATIO-readslab-${var} "" test_readslab_${var}.out test_mat readslab test_mat_writeslab.mat ${var})
#ENDFOREACH()

#FOREACH(arg write_struct_2d_numeric write_struct_complex_2d_numeric)
#    SET(VERSIONS 4 5)
#    IF (MAT73)
#        SET(VERSIONS ${VERSIONS} 7.3)
#    ENDIF()
#    FOREACH(field 1 2)
#        FOREACH(class double single int64 uint64 int32 uint32 int16 uint16 int8 uint8)
#            FOREACH (version ${VERSIONS})
#                SET(testname ${arg}-${class}-${version})
#                SET(filename test_${testname}.mat)
#                SET(testname getstructfield-${testname})
#                MATIO_TEST(MATIO-${testname} "" ${testname}.out test_mat getstructfield ${filename} a field${field})
#                IF (${version} STREQUAL "5")
#                    SET(testname ${arg}-${class}-${version}-compressed)
#                    SET(filename test_${testname}.mat)
#                    SET(testname getstructfield-${testname})
#                    MATIO_TEST(MATIO-${testname} "" ${testname}.out test_mat getstructfield ${filename} a field${field})
#                ENDIF() 
#            ENDFOREACH()
#        ENDFOREACH()
#    ENDFOREACH()
#ENDFOREACH()

# Add more tests for these.

SET(MATIO_READ_TESTS readvar4 readslab4 slab3)

# Set tests that are expected to fail (TO BE CORRECTED).

IF (ENABLE_FORTRAN)
    INCLUDE_DIRECTORIES(${MATIO_SOURCE_DIR}/src/fortran ${MATIO_BINARY_DIR}/src/fortran)
    ADD_EXECUTABLE(test_matf test_matf.f90)
    TARGET_LINK_LIBRARIES(test_matf matio m)
    # TESTS
ENDIF()
