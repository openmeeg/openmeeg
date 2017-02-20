#  Adapted from Caffe
find_package(OpenBLAS QUIET CONFIG)

if (OpenBLAS_FOUND) # the git version propose a OpenBLASConfig.cmake
    message(STATUS "OpenBLASConfig found")
    set(OpenBLAS_INCLUDE_DIR ${OpenBLAS_INCLUDE_DIRS})
    set(OpenBLAS_LIB ${OpenBLAS_LIBRARIES})
else()
    SET(OpenBLAS_INCLUDE_SEARCH_PATHS
        /usr/include
        /usr/include/openblas
        /usr/include/openblas-base
        /usr/local/include
        /usr/local/include/openblas
        /usr/local/include/openblas-base
        /opt/OpenBLAS/include
        $ENV{OpenBLAS_HOME}
        $ENV{OpenBLAS_HOME}/include
        )

    SET(OpenBLAS_LIB_SEARCH_PATHS
        /lib/
        /lib/openblas-base
        /lib64/
        /usr/lib
        /usr/lib/openblas-base
        /usr/lib64
        /usr/local/lib
        /usr/local/lib64
        /opt/OpenBLAS/lib
        $ENV{OpenBLAS}cd
        $ENV{OpenBLAS}/lib
        $ENV{OpenBLAS_HOME}
        $ENV{OpenBLAS_HOME}/lib
        )

    FIND_PATH(OpenBLAS_INCLUDE_DIR NAMES openblas_config.h PATHS ${OpenBLAS_INCLUDE_SEARCH_PATHS})
    FIND_LIBRARY(OpenBLAS_LIB NAMES openblas PATHS ${OpenBLAS_LIB_SEARCH_PATHS})
    # mostly for debian
    FIND_LIBRARY(Lapacke_LIB NAMES lapacke PATHS ${OpenBLAS_LIB_SEARCH_PATHS})

    SET(OpenBLAS_FOUND ON)

    #    Check include files
    IF(NOT OpenBLAS_INCLUDE_DIR)
        SET(OpenBLAS_FOUND OFF)
        MESSAGE(STATUS "Could not find OpenBLAS include. Turning OpenBLAS_FOUND off")
    ENDIF()

    #    Check libraries
    IF(NOT OpenBLAS_LIB)
        SET(OpenBLAS_FOUND OFF)
        MESSAGE(STATUS "Could not find OpenBLAS lib. Turning OpenBLAS_FOUND off")
    ENDIF()

    IF (OpenBLAS_FOUND)
        SET(OpenBLAS_LIBRARIES ${OpenBLAS_LIB})
        IF (NOT Lapacke_LIB-NOTFOUND)
            SET(OpenBLAS_LIBRARIES ${OpenBLAS_LIBRARIES} ${Lapacke_LIB})
        endif()
        IF (NOT OpenBLAS_FIND_QUIETLY)
            MESSAGE(STATUS "Found OpenBLAS libraries: ${OpenBLAS_LIBRARIES}")
            MESSAGE(STATUS "Found OpenBLAS include: ${OpenBLAS_INCLUDE_DIR}")
        ENDIF (NOT OpenBLAS_FIND_QUIETLY)
    ELSE (OpenBLAS_FOUND)
        IF (OpenBLAS_FIND_REQUIRED)
            MESSAGE(FATAL_ERROR "Could not find OpenBLAS")
        ENDIF (OpenBLAS_FIND_REQUIRED)
    ENDIF (OpenBLAS_FOUND)

endif()

MARK_AS_ADVANCED(
    OpenBLAS_INCLUDE_DIR
    OpenBLAS_LIB
    OpenBLAS
    )
