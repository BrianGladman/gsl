noinst_LTLIBRARIES = libgslmultimin.la 

pkginclude_HEADERS = gsl_multimin.h

noinst_HEADERS = 

INCLUDES= -I$(top_builddir)

libgslmultimin_la_SOURCES = directional_minimize.c fdfminimizer.c steepest_descent.c conjugate.c convergence.c vector_bfgs.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_funcs.c test_funcs.h
test_LDADD = libgslmultimin.la ../min/libgslmin.la ../blas/libgslblas.la ../blas/libgslblasnative.la ../blas/libgslblascblas.la ../linalg/libgsllinalg.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la
test_DEPENDENCIES = libgslmultimin.la

