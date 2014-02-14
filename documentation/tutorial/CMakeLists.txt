INCLUDE(UseLATEX.cmake)

ADD_LATEX_DOCUMENT(tutorial.tex
                   INPUTS tutorial.tex
                   IMAGE_DIRS images
                   CONFIGURE tutorial.tex
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES ON)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/tutorial.pdf DESTINATION share/doc/OpenMEEG/tutorials)
