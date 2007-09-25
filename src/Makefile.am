# MatLib settings

INCLUDES= -I$(top_srcdir) -I$(top_srcdir)/src \
          -I$(top_srcdir)/libs/MatLib/MatLib $(BLAS_INCLUDES) 

LDADD= $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

AM_CXXFLAGS= -O3 -finline-functions -ffast-math -fstrength-reduce -frerun-loop-opt \
             -fgcse -fexpensive-optimizations -funroll-loops -Winline

bin_PROGRAMS= forward gain assemble inverse matbin2txt mattxt2bin mesh_convert \
              minverser presmooth symbin2txt symtxt2bin txt2tex

forward_SOURCES= forward.cpp
gain_SOURCES= gain.cpp
assemble_SOURCES= assembleFerguson.cpp assembleLHS.cpp assembleRHS.cpp assembleSensors.cpp \
                  danielsson.cpp geometry.cpp assemble.cpp mesh3.cpp operators.cpp
inverse_SOURCES= inverse.cpp
matbin2txt_SOURCES= matbin2txt.cpp
mattxt2bin_SOURCES= mattxt2bin.cpp
mesh_convert_SOURCES= mesh_convert.cpp
minverser_SOURCES= minverser.cpp
presmooth_SOURCES= presmooth.cpp
symbin2txt_SOURCES= symbin2txt.cpp
symtxt2bin_SOURCES= symtxt2bin.cpp
txt2tex_SOURCES= txt2tex.cpp
