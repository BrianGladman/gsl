noinst_LIBRARIES = libroots.a 

noinst_PROGRAMS = test

include_HEADERS = gsl_roots.h

INCLUDES = -I$(srcdir)/../err

TESTS = test

test_LDADD = libroots.a ../err/liberr.a

test_DEPENDENCIES = libroots.a

test_SOURCES = test.c

libroots_a_SOURCES = bisection.c
