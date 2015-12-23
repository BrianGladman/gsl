noinst_LTLIBRARIES = libgslmultifit_nlinear.la 

pkginclude_HEADERS = gsl_multifit_nlinear.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslmultifit_nlinear_la_SOURCES = fdf.c

noinst_HEADERS =

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslmultifit_nlinear.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../sort/libgslsort.la ../statistics/libgslstatistics.la ../vector/libgslvector.la ../block/libgslblock.la  ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la
