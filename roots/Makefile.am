# -*-makefile-*-

noinst_LTLIBRARIES = libgslroots.la 

pkginclude_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES= -I$(top_builddir)

libgslroots_la_SOURCES = bisection.c brent.c falsepos.c \
			newton.c  secant.c steffenson.c \
			convergence.c fsolver.c fdfsolver.c 

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_roots.c test_funcs.c test_macros.c test.h
test_LDADD = libgslroots.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la
test_DEPENDENCIES = libgslroots.la
