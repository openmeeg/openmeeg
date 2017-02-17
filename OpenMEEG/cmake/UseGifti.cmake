#------------------------------------------------------------
# GIFTI C library
#------------------------------------------------------------

option(USE_GIFTI "Build the project using GIFTI IO support" OFF)
mark_as_advanced(USE_GIFTI)

if (USE_GIFTI)
    find_package(EXPAT)
    find_package(ZLIB)
    find_library(NIFTI_LIBRARY niftiio)
    find_library(GIFTI_LIBRARY giftiio)
    find_library(ZNZ_LIBRARY znz)
    set(NIFTI_LIBRARIES ${EXPAT_LIBRARIES} ${ZLIB_LIBRARIES} ${NIFTI_LIBRARY} ${ZNZ_LIBRARY} m)
    find_path(GIFTI_INCLUDE_PATH gifti_io.h PATHS /usr/include/gifti)
    include_directories(${GIFTI_INCLUDE_PATH})
    find_path(NIFTI_INCLUDE_PATH nifti1_io.h PATHS /usr/include/nifti)
    include_directories(${NIFTI_INCLUDE_PATH})
    set(OPENMEEG_OTHER_INCLUDE_DIRECTORIES ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES} ${GIFTI_INCLUDE_PATH} ${NIFTI_INCLUDE_PATH})
    set(OPENMEEG_LIBRARIES ${OPENMEEG_LIBRARIES} ${GIFTI_LIBRARY} ${NIFTI_LIBRARIES})
endif()
