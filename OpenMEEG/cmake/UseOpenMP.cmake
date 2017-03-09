#------------------------------------------------------------
# OpenMP library
#------------------------------------------------------------

option(USE_OMP "Use OpenMP" OFF)

if (USE_OMP)
    find_package(OpenMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DUSE_OMP)
endif()
