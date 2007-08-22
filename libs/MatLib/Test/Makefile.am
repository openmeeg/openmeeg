INCLUDES= -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/libs/MatLib/MatLib $(BLAS_INCLUDES)
LDADD= $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

noinst_PROGRAMS= test
test_SOURCES= test.cpp chrono.cpp
