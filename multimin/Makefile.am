pkglib_LIBRARIES = libgslmultimin.a 

pkginclude_HEADERS = gsl_multimin.h

noinst_HEADERS = 

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../linalg -I$(srcdir)/../matrix -I$(srcdir)/../vector -I$(srcdir)/../block -I$(srcdir)/../blas -I$(srcdir)/../err -I$(srcdir)/../min

libgslmultimin_a_SOURCES = directional_minimize.c fdfminimizer.c steepest_descent.c conjugate.c convergence.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c 
test_LDADD = libgslmultimin.a ../min/libgslmin.a ../blas/libgslblas.a ../blas/libgslblasnative.a ../blas/libgslblascblas.a ../linalg/libgsllinalg.a ../matrix/libgslmatrix.a ../vector/libgslvector.a ../block/libgslblock.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslmultimin.a

