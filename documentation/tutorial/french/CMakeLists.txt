ADD_LATEX_DOCUMENT(tutoriel.tex
                       INPUTS tutoriel.tex Contexte.tex Commandes.tex Annexes.tex Donnees.tex
                       IMAGE_DIRS ../images
                       CONFIGURE tutoriel.tex
                       DEFAULT_PDF MANGLE_TARGET_NAMES]
)