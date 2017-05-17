# install headers

set(OpenMEEGMaths_HEADERS
    fast_sparse_matrix.h linop.h OpenMEEGMathsConfig.h 
    matrix.h symmatrix.h sparse_matrix.h vector.h OMassert.H
    #   These files are imported from another repository.
    #   Please do not update them in this repository.
    AsciiIO.H BrainVisaTextureIO.H Exceptions.H IOUtils.H MathsIO.H MatlabIO.H RC.H 
    TrivialBinIO.H)

install(FILES ${OpenMEEGMaths_HEADERS}
    DESTINATION ${OpenMEEG_HEADER_INSTALLDIR} COMPONENT Development)

install(DIRECTORY BlasLapackImplementations
    DESTINATION ${OpenMEEG_HEADER_INSTALLDIR} COMPONENT Development)
