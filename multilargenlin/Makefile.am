noinst_LTLIBRARIES = libgslmultilargenlin.la 

pkginclude_HEADERS = gsl_multilarge_nlin.h

libgslmultilargenlin_la_SOURCES = lmnielsen.c nlinear.c

AM_CPPFLAGS = -I$(top_srcdir)

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmultilargenlin.la ../test/libgsltest.la ../multilarge/libgslmultilarge.la ../multifit/libgslmultifit.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../sys/libgslsys.la ../utils/libutils.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../complex/libgslcomplex.la
