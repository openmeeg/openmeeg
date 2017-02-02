set(Atlas    USE_ATLAS)
set(OpenBlas USE_OPENBLAS)
set(MKL      USE_MKL)
set(Lapack   USE_LAPACK)
set(ACML     USE_ACML)

if (BLASLAPACK_IMPLEMENTATION)
    set(BLASLAPACK_IMPLEMENTATION_DEFAULT ${BLASLAPACK_IMPLEMENTATION})
else()
    if (WIN32)
        set(BLASLAPACK_IMPLEMENTATION_DEFAULT MKL)
    else()
        set(BLASLAPACK_IMPLEMENTATION_DEFAULT Lapack)
    endif()
endif()

set (BLASLAPACK_IMPLEMENTATION ${BLASLAPACK_IMPLEMENTATION_DEFAULT} CACHE STRING "Choose 
                the proper Blas/Lapack implementation Atlas/OpenBlas/MKL/ACML/Lapack" FORCE)
# Set the possible values of build type for cmake-gui
set_property(CACHE BLASLAPACK_IMPLEMENTATION PROPERTY STRINGS "Atlas" "OpenBlas" "MKL" "ACML" "Lapack")

if (NOT ${BLASLAPACK_IMPLEMENTATION})
    message(ERROR "Unknown blas/lapack implementation in BLASLAPACK_IMPLEMENTATION")
endif()

#   Ensure that only one lapack implementation is selected by clearing all variable before setting the one chosen.

foreach (i Atlas OpenBlas MKL Lapack)
    unset(${${i}})
endforeach()
set(${${BLASLAPACK_IMPLEMENTATION}} ON)
