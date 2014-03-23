# install headers

set(OPENMEEGMATHS_HEADERS 
    DLLDefinesOpenMEEGMaths.h fast_sparse_matrix.h linop.h OpenMEEGMathsConfig.h 
    matrix.h symmatrix.h sparse_matrix.h vector.h
    #   These files are imported from another repository.
    #   Please do not update them in this repository.
    AsciiIO.H BrainVisaTextureIO.H Exceptions.H IOUtils.H MathsIO.H MatlabIO.H RC.H 
    TrivialBinIO.H)

install(FILES ${OPENMEEGMATHS_HEADERS}
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)
