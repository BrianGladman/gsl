noinst_LTLIBRARIES = libgslmultilarge.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_multilarge.h

libgslmultilarge_la_SOURCES = multilarge.c normal.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslmultilarge.la ../test/libgsltest.la ../multifit/libgslmultifit.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../utils/libutils.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la

test_SOURCES = test.c
