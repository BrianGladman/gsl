pkglib_LIBRARIES = libgslmin.a 

pkginclude_HEADERS = gsl_min.h

noinst_HEADERS = min.h

INCLUDES= -I$(top_builddir)

libgslmin_a_SOURCES = fsolver.c golden.c brent.c convergence.c bracketing.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_funcs.c test.h
test_LDADD = libgslmin.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a ../sys/libgslsys.a
test_DEPENDENCIES = libgslmin.a
