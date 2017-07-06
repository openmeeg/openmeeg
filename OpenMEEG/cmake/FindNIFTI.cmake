
find_package(EXPAT ${FIND_MODE})
find_package(ZLIB ${FIND_MODE})
find_library(NIFTI_LIBRARY niftiio ${FIND_MODE})

find_library(ZNZ_LIBRARY znz ${FIND_MODE})
find_path(NIFTI_INCLUDE_PATH nifti1_io.h PATH_SUFFIXES nifti ${FIND_MODE})

if (NIFTI_FIND_REQUIRED)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(NIFTI DEFAULT_MSG NIFTI_INCLUDE_PATH NIFTI_LIBRARY)
endif()
