#------------------------------------------------------------
# Atlas library
#------------------------------------------------------------

if (USE_ATLAS)
    find_package(Atlas ${FIND_MODE} MODULE)
    if (Atlas_FOUND)
        include_directories(${Atlas_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${Atlas_LIBRARIES})
        list(APPEND OpenMEEG_DEPENDENCIES Atlas)
    endif()
endif()
