pkglib_LIBRARIES = libgslmultimin.a 

pkginclude_HEADERS = gsl_multimin.h

noinst_HEADERS = 

INCLUDES= -I$(top_builddir)

libgslmultimin_a_SOURCES = directional_minimize.c fdfminimizer.c steepest_descent.c conjugate.c convergence.c vector_bfgs.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_funcs.c test_funcs.h
test_LDADD = libgslmultimin.a ../min/libgslmin.a ../blas/libgslblas.a ../blas/libgslblasnative.a ../blas/libgslblascblas.a ../linalg/libgsllinalg.a ../matrix/libgslmatrix.a ../vector/libgslvector.a ../block/libgslblock.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslmultimin.a

