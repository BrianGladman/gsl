noinst_LTLIBRARIES = libgslfit.la 

pkginclude_HEADERS = gsl_fit.h

INCLUDES= -I$(top_builddir)

libgslfit_la_SOURCES = linear.c multiple.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c
test_LDADD = libgslfit.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la
test_DEPENDENCIES = libgslfit.la
