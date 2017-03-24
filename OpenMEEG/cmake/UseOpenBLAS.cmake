#------------------------------------------------------------
# OpenBLAS library
#------------------------------------------------------------

if (USE_OPENBLAS)
    if(WIN32)
        # download Openblas
        include(ExternalProject)
        ExternalProject_Add(OpenBLAS DOWNLOAD_NAME OpenBLAS URL "https://downloads.sourceforge.net/project/openblas/v0.2.19/OpenBLAS-v0.2.19-Win64-int32.zip" URL_MD5 "7ff6092397a93494c137e23670dd72ec")
        set(OpenBLAS_INCLUDE_DIR OpenBLAS/include)
        set(OpenBLAS_LIBRARIES OpenBLAS/lib/libopenblas.a)
        message("------------------------------------------------------------------------------ ${OpenBLAS_INCLUDE_DIR} ")
        message("------------------------------------------------------------------------------ ${OpenBLAS_LIBRARIES} ")
        set(OpenBLAS_FOUND)
    else()
        find_package(OpenBLAS ${REQUIRED} MODULE)
    endif()
    if (OpenBLAS_FOUND)
        include_directories(${OpenBLAS_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
    endif()
endif()
