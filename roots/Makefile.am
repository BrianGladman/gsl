# -*-makefile-*-
# $Id$

noinst_LIBRARIES = libroots.a 

noinst_PROGRAMS = test

include_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES = -I$(srcdir)/../err

TESTS = test-macros test-bisection test-falsepos test-secant test-newton

test_LDADD = libroots.a ../err/liberr.a

test_DEPENDENCIES = libroots.a

test_SOURCES = test.c test.h

libroots_a_SOURCES = utility.c bisection.c falsepos.c secant.c newton.c
