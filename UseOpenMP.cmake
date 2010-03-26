OPTION(USE_OMP "Use OpenMP" OFF)

IF(USE_OMP)
    IF(UNIX)
        IF(USE_ICC)
            #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -axT")
            #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp -axT")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")
        ELSE(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
            #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
            #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse")
        ENDIF(USE_ICC)
    ELSE(UNIX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
    ENDIF(UNIX)
    ADD_DEFINITIONS(-DUSE_OMP)
ENDIF(USE_OMP)
