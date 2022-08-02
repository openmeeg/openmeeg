set(CHECK_CXX_FEATURE_PREFIX "")
mark_as_advanced(CHECK_CXX_FEATURE_PREFIX)

macro(CHECK_CXX_FEATURE feature file message)
    if (NOT ${feature}_DONE)
        set(flags)
        foreach (flag ${ARGV3})
            set(flags ${flags} -D${flag})
        endforeach()
        set(libs ${ARGV4})
        message(STATUS "Check whether the compiler ${message}")
        try_compile(RESULT ${CMAKE_BINARY_DIR}
            #${CMAKE_ROOT}/Modules/TestForSTDNamespace.cxx
            ${CMAKE_SOURCE_DIR}/cmake/utils/cxx_tests/${file}
            COMPILE_DEFINITIONS "${CHECK_CXX_FEATURE_DEFINITIONS} ${flags}"
            LINK_LIBRARIES "${CHECK_CXX_FEATURE_LINK_LIBRARIES} ${libs}"
            OUTPUT_VARIABLE OUTPUT)

        if (RESULT)
            set (FOUND "found")
            set (STATUS "passed")
            set(CHECK_CXX_FEATURE_DEFINITIONS "${CHECK_CXX_FEATURE_DEFINITIONS} -D${feature}")
        else()
            set (FOUND "not found")
            set (STATUS "failed")
        endif()
        message(STATUS "Check whether the compiler ${message} - ${FOUND}")
        set("${CHECK_CXX_FEATURE_PREFIX}${feature}" ${RESULT} CACHE INTERNAL "Does the compiler ${message}")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
             "Determining if the CXX compiler ${message} ${STATUS} with "
             "the following output:\n${OUTPUT}\n\n")
    endif()
    set(${feature}_DONE 1)
endmacro()

macro(CHECK_CXX_OPENMP_SUPPORT)
    CHECK_CXX_FEATURE(HAVE_OPENMP_RANGEFOR openmp_support.cpp  "has OpenMP support for range for loops" RANGEFOR ${OpenMP_CXX_LIBRARIES})
    CHECK_CXX_FEATURE(HAVE_OPENMP_ITERATOR openmp_support.cpp  "has OpenMP support for iterator loops"  ITERATOR ${OpenMP_CXX_LIBRARIES})
    CHECK_CXX_FEATURE(HAVE_OPENMP_UNSIGNED openmp_support.cpp  "has OpenMP support for unsigned loops"  ""       ${OpenMP_CXX_LIBRARIES})
endmacro()
