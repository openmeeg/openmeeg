#------------------------------------------------------------
# GIFTI C library
#------------------------------------------------------------

option(USE_GIFTI "Use GIFTI IO support" OFF)

if (USE_GIFTI)
    find_package(EXPAT)
    find_package(ZLIB)
    find_library(NIFTI_LIBRARY niftiio)
    find_library(GIFTI_LIBRARY giftiio)
    find_library(ZNZ_LIBRARY znz)
    find_path(GIFTI_INCLUDE_PATH gifti_io.h PATH_SUFFIXES gifti)
    find_path(NIFTI_INCLUDE_PATH nifti1_io.h PATH_SUFFIXES nifti)
    set(GIFTI_INCLUDE_DIRS ${GIFTI_INCLUDE_PATH} ${NIFTI_INCLUDE_PATH})
    set(GIFTI_LIBRARIES ${EXPAT_LIBRARIES} ${ZLIB_LIBRARIES} ${NIFTI_LIBRARY} ${ZNZ_LIBRARY} m ${NIFTI_LIBRARY} ${GIFTI_LIBRARY})
    list(APPEND OpenMEEG_OTHER_INCLUDE_DIRS ${GIFTI_INCLUDE_DIRS})
    set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${GIFTI_LIBRARIES})
endif()
