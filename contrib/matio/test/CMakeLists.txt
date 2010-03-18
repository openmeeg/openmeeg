set(TESTS test_mat test_snprintf)

FOREACH(test ${TESTS})
    ADD_EXECUTABLE(${test} ${test}.c)
    TARGET_LINK_LIBRARIES(${test} matio m)
ENDFOREACH()

ADD_TEST(test_snprintf test_snprintf)

SET(test_mat_args
    copy readvar getstructfield readvarinfo readslab write writecompressed writesparse write_compressed_sparse
    write_struct write_compressed_struct writecell write_compressed_cell writeinf writenan writenull writeslab
    readvar4 ind2sub sub2ind)

FOREACH(arg ${test_mat_args})
    ADD_TEST(test_mat_${arg} test_mat ${arg})
ENDFOREACH()

IF (ENABLE_FORTRAN)
    INCLUDE_DIRECTORIES(${MATIO_SOURCE_DIR}/src/fortran ${MATIO_BINARY_DIR}/src/fortran)
    ADD_EXECUTABLE(test_matf test_matf.f90)
    TARGET_LINK_LIBRARIES(test_matf matio m)
    # TESTS
ENDIF()
