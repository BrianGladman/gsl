# -*-makefile-*-

noinst_LIBRARIES = libgslroots.a 

check_PROGRAMS = test

include_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES = -I$(srcdir)/../err

TESTS = test

test_SOURCES = test.c test.h
test_LDADD = libgslroots.a ../err/libgslerr.a ../err/libgsltest.a
test_DEPENDENCIES = libgslroots.a

libgslroots_a_SOURCES = bisection.c falsepos.c secant.c newton.c
