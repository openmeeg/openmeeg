ADD_LATEX_DOCUMENT(tutorial.tex
                   INPUTS tutorial.tex Context.tex Commands.tex Appendix.tex Data.tex
                   IMAGE_DIRS ../images
                   CONFIGURE tutorial.tex
                   DEFAULT_PDF MANGLE_TARGET_NAMES])

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/tutorial.pdf DESTINATION share/doc/OpenMEEG/tutorials)
