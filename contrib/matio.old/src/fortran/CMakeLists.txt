INCLUDE(FortranCInterface)
FortranCInterface_HEADER(fortran.h MACRO_NAMESPACE "FC_" FROM_LIBRARY blas[daxpy])

SET(SIZEOF_VOID_P ${SIZEOF_VOID_PTR})
Configure_File(matio_t.inc.in matio_t.inc IMMEDIATE @ONLY)
