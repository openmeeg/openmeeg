#---------------------------------------------------------------
# Test and setup the C++ compiler features
#---------------------------------------------------------------

option(USE_OPENMP "Use OpenMP (if available)"  ON)

if(CMAKE_CXX_COMPILER_LOADED)
     message(STATUS "Compiler path: ${CMAKE_CXX_COMPILER}")
     message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
     message(STATUS "Compiler version:
             ${CMAKE_CXX_COMPILER_VERSION}")
endif()

if (USE_OPENMP)
    include(CheckCXXFeatures)
    if (OPENMP_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES_OLD ${CMAKE_FIND_LIBRARY_SUFFIXES})
        if (NOT WIN32)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        else()
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
        endif()
    endif()
    find_package(OpenMP)
    if (OPENMP_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_OLD})
    endif()
    if (OpenMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        CHECK_CXX_OPENMP_SUPPORT()
        set(OPENMP_DEFINITIONS USE_OMP)
        foreach (i OPENMP_RANGEFOR OPENMP_ITERATOR OPENMP_UNSIGNED)
            if (HAVE_${i})
                set(OPENMP_DEFINITIONS ${OPENMP_DEFINITIONS} ${i})
            endif()
        endforeach()
        message("-- Found OpenMP library " ${OpenMP_CXX_LIBRARIES})
    endif()
endif()

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "None" "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

if (WIN32 AND NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
endif()
