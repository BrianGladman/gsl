# -*-makefile-*-

pkglib_LIBRARIES = libgslroots.a 

pkginclude_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES= -I$(top_builddir)

libgslroots_a_SOURCES = bisection.c brent.c falsepos.c \
			newton.c  secant.c steffenson.c \
			convergence.c fsolver.c fdfsolver.c 

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_roots.c test_funcs.c test_macros.c test.h
test_LDADD = libgslroots.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslroots.a
