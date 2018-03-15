
if (USE_AUTO)
    # this macro tries to use library in the defined priority
    macro(priorities)
        foreach (lib ${ARGN})
            string(TOUPPER ${lib} LIB)
            set(USE_${LIB} true)
            include(Use${lib})
            if (${lib}_FOUND)
                # set the found BLASLAPACK_IMPLEMENTATION (in case it was Auto)
                set(BLASLAPACK_IMPLEMENTATION "${lib}" CACHE STRING "${BLASLAPACK_IMPLEMENTATION_DOCSTRING}" FORCE)
                break()
            else()
                unset(USE_${LIB})
            endif()
        endforeach()
    endmacro()

    set(FIND_MODE "QUIET")

    # define here the priorities of BLAS/LAPACK to use
    if (WIN32)
        priorities(MKL OpenBLAS LAPACK)
    elseif (APPLE)
        priorities(MKL OpenBLAS vecLib LAPACK)
    else()
        priorities(MKL OpenBLAS Atlas LAPACK)
    endif()
endif()
