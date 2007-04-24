# MatLib settings

INCLUDES= -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/libs/LsLib \
          -I$(top_srcdir)/libs/MatLib/MatLib $(BLAS_INCLUDES) 

LDADD= $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

AM_CXXFLAGS= -O3 -finline-functions -ffast-math -fstrength-reduce -frerun-loop-opt \
             -fgcse -fexpensive-optimizations -funroll-loops -Winline

bin_PROGRAMS= forward gain assemble inverse matbin2txt mattxt2bin meshconvert \
              minverser presmooth symbin2txt symtxt2bin txt2tex

forward_SOURCES= forward.cpp
gain_SOURCES= gain.cpp
assemble_SOURCES= assembleFerguson.cpp assembleMat.cpp assemble_RHSmatrix.cpp assemble_squids_EEG.cpp \
                  danielsson.cpp geometry.cpp assemble.cpp mesh3.cpp operateurs.cpp
inverse_SOURCES= inverse.cpp
matbin2txt_SOURCES= matbin2txt.cpp
mattxt2bin_SOURCES= mattxt2bin.cpp
meshconvert_SOURCES= meshconvert.cpp
minverser_SOURCES= minverser.cpp
presmooth_SOURCES= presmooth.cpp
symbin2txt_SOURCES= symbin2txt.cpp
symtxt2bin_SOURCES= symtxt2bin.cpp
txt2tex_SOURCES= txt2tex.cpp
