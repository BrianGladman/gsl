noinst_LTLIBRARIES = libgslmultifit.la 

pkginclude_HEADERS = gsl_multifit.h

INCLUDES= -I$(top_builddir)

libgslmultifit_la_SOURCES = multilinear.c 

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c
test_LDADD = libgslmultifit.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la
test_DEPENDENCIES = libgslmultifit.la
