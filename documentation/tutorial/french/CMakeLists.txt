#ADD_LATEX_DOCUMENT(tutorial-french.tex
#                   INPUTS tutorial-french.tex Contexte.tex Commandes.tex Annexes.tex Donnees.tex
#                   IMAGE_DIRS ../images
#                   CONFIGURE tutorial-french.tex
#                   DEFAULT_PDF
#                   MANGLE_TARGET_NAMES ON)

#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/tutorial-french.pdf DESTINATION share/doc/OpenMEEG/tutorials)
