# -*-makefile-*-

pkglib_LIBRARIES = libgslmultiroots.a 

pkginclude_HEADERS = gsl_multiroots.h

noinst_HEADERS = enorm.c dogleg.c

INCLUDES= -I$(top_builddir)

libgslmultiroots_a_SOURCES = fdjac.c fsolver.c fdfsolver.c convergence.c  newton.c gnewton.c dnewton.c broyden.c hybrid.c hybridj.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_funcs.c test_funcs.h
test_LDADD = libgslmultiroots.a ../linalg/libgsllinalg.a ../matrix/libgslmatrix.a ../vector/libgslvector.a ../block/libgslblock.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslmultiroots.a
