noinst_LTLIBRARIES = libgslsplinalg.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_splinalg.h

libgslsplinalg_la_SOURCES = spgmres.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslsplinalg.la ../spblas/libgslspblas.la ../test/libgsltest.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../utils/libutils.la ../rng/libgslrng.la ../err/libgslerr.la

test_SOURCES = test.c
