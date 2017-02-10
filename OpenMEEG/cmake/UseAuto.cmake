
if (USE_AUTO)
    # this macro tries to use library in the defined priority
    macro(priorities)
        foreach (lib ${ARGN})
            string(TOUPPER ${lib} LIB)
            set(USE_${LIB} true)
            include(Use${lib})
            if (${lib}_FOUND)
                break()
            else()
                unset(USE_${LIB})
            endif()
        endforeach()
    endmacro()

    set(REQUIRED "QUIET")

    # define here the priorities of BLAS/LAPACK to use
    if (WIN32)
        priorities(MKL OpenBLAS LAPACK)
    elseif (APPLE)
        priorities(OpenBLAS vecLib LAPACK)
    else()
        priorities(OpenBLAS Atlas LAPACK)
    endif()
endif()
