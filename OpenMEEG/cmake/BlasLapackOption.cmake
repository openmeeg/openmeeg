
set(Atlas    USE_ATLAS)
set(OpenBLAS USE_OPENBLAS)
set(MKL      USE_MKL)
set(LAPACK   USE_LAPACK)
set(vecLib   USE_VECLIB)
set(Auto     USE_AUTO)

# the default case
set(BLASLAPACK_IMPLEMENTATION_DEFAULT Auto)
if (BLASLAPACK_IMPLEMENTATION)
    set(BLASLAPACK_IMPLEMENTATION_DEFAULT ${BLASLAPACK_IMPLEMENTATION})
endif()

# the list of possibilites depending on the OS
set(LIST_IMPL "Auto" "MKL" "LAPACK" "OpenBLAS")
if (APPLE)
    set(LIST_IMPL ${LIST_IMPL} "vecLib")
else()
    set(LIST_IMPL ${LIST_IMPL} "Atlas")
endif()

set (BLASLAPACK_IMPLEMENTATION ${BLASLAPACK_IMPLEMENTATION_DEFAULT} CACHE STRING 
    "Choose the proper Blas/Lapack implementation: ${LIST_IMPL}" FORCE)

# Set the possible values of build type for cmake-gui
set_property(CACHE BLASLAPACK_IMPLEMENTATION PROPERTY STRINGS ${LIST_IMPL})

# Ensure that only one lapack implementation is selected by clearing all variable before setting the one chosen.
foreach (i Auto Atlas OpenBLAS MKL LAPACK vecLib)
    unset(${${i}})
endforeach()
set(${${BLASLAPACK_IMPLEMENTATION}} ON)

#unset unused variables
foreach (i Auto Atlas OpenBLAS MKL LAPACK vecLib)
    unset(${i})
endforeach()
