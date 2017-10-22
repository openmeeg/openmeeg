#------------------------------------------------------------
# GIFTI C library
#------------------------------------------------------------

option(USE_GIFTI "Use GIFTI IO support" OFF)

if (USE_GIFTI)
    set(FIND_MODE "QUIET")
    if(CMAKE_PROJECT_NAME STREQUAL "OpenMEEG")
        set(FIND_MODE "REQUIRED")
    endif()

    find_package(NIFTI ${FIND_MODE} MODULE)
    find_package(GIFTI ${FIND_MODE} MODULE)

    set(GIFTI_INCLUDE_DIRS ${GIFTI_INCLUDE_PATH} ${NIFTI_INCLUDE_PATH})
    set(GIFTI_LIBRARIES ${EXPAT_LIBRARIES} ${ZLIB_LIBRARIES} ${NIFTI_LIBRARY} ${ZNZ_LIBRARY} m ${NIFTI_LIBRARY} ${GIFTI_LIBRARY})
    list(APPEND OpenMEEG_OTHER_INCLUDE_DIRS ${GIFTI_INCLUDE_DIRS})
    set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${GIFTI_LIBRARIES})
endif()
