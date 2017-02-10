
if (USE_ATLAS)
    find_package(Atlas ${REQUIRED} MODULE)
    if (Atlas_FOUND)
        include_directories(${Atlas_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${Atlas_LIBRARIES})
    endif()
endif()
