noinst_LTLIBRARIES = libgslinterpolation.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_interp.h

noinst_HEADERS = bsearch.h


INCLUDES= -I$(top_builddir)

TESTS = test

test_LDADD = libgslinterpolation.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../blas/libgslblasnative.la ../vector/libgslvector.la ../block/libgslblock.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

test_DEPENDENCIES = libgslinterpolation.la
test_SOURCES = test.c

libgslinterpolation_la_SOURCES = accel.c akima.c bsearch.c cspline.c interp.c linear.c integ_eval_macro.h
