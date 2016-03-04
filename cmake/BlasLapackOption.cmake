set(Atlas    USE_ATLAS)
set(OpenBlas USE_OPENBLAS)
set(MKL      USE_MKL)
set(Lapack   USE_LAPACK)
set(ACML     USE_ACML)

if (WIN32)
    set(BLASLAPACK_IMPLEMENTATION_DEFAULT MKL)
else()
    set(BLASLAPACK_IMPLEMENTATION_DEFAULT Atlas)
endif()

set(BLASLAPACK_IMPLEMENTATION "OpenBlas" CACHE FILEPATH "Vhoose the proper Blas/Lapack implementation Atlas/OpenBlas/MKL/ACML/Lapack")
if (NOT ${BLASLAPACK_IMPLEMENTATION})
    message(ERROR "Unknown blas/lapack implementation in BLASLAPACK_IMPLEMENTATION")
endif()

#   Ensure that only one lapack implementation is selected by clearing all variable before setting the one chosen.

foreach (i Atlas OpenBlas MKL Lapack)
    unset(${${i}})
endforeach()
set(${${BLASLAPACK_IMPLEMENTATION}} ON)
