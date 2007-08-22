# MatLib settings

INCLUDES= -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/libs/MatLib/MatLib $(BLAS_INCLUDES)
LDADD= $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

AM_CXXFLAGS= -O3 -finline-functions -ffast-math -fstrength-reduce -frerun-loop-opt -fgcse -fexpensive-optimizations -funroll-loops -Winline

bin_PROGRAMS= load_geo inverse_test
load_geo_SOURCES= load_geo.cpp
inverse_test_SOURCES= inverse_test.cpp
