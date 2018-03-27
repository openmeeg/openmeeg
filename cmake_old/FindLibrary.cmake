
option(BUILD_SHARED_LIBS "Build shared libs" ON)
mark_as_advanced(BUILD_SHARED_LIBS)

if (APPLE_STANDALONE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so;.dylib")
endif()

string(COMPARE NOTEQUAL "${BUILD_SHARED_STATUS}" "" BUILD_SHARED_STATUS_NOT_EMPTY)
if (BUILD_SHARED_STATUS_NOT_EMPTY)
    string(COMPARE NOTEQUAL "${BUILD_SHARED_STATUS}" "${BUILD_SHARED_LIBS}" RESET)
endif()

# Store in cache previous value of BUILD_SHARED_LIBS
set(BUILD_SHARED_STATUS "${BUILD_SHARED_LIBS}" CACHE INTERNAL "Previous shared status" FORCE)

function(FIND_LIBRARY VAR)
    if (${RESET})
        set(${VAR} NOTFOUND CACHE STRING "" FORCE)
    endif()
    _FIND_LIBRARY(${VAR} ${ARGN})
    mark_as_advanced(${VAR})
endfunction()
