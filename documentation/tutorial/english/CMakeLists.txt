ADD_LATEX_DOCUMENT(tutorial-english.tex
                   INPUTS tutorial-english.tex Context.tex Commands.tex Appendix.tex Data.tex
                   IMAGE_DIRS ../images
                   CONFIGURE tutorial-english.tex
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES ON)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/tutorial-english.pdf DESTINATION share/doc/OpenMEEG/tutorials)
