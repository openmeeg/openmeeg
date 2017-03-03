#------------------------------------------------------------
# OpenBLAS library
#------------------------------------------------------------

if (USE_OPENBLAS)
    find_package(OpenBLAS ${REQUIRED} MODULE)
    if (OpenBLAS_FOUND)
        include_directories(${OpenBLAS_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
    endif()
endif()
