# -*-makefile-*-

lib_LIBRARIES = libgslroots.a 

noinst_PROGRAMS = test

include_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES = -I$(srcdir)/../err

TESTS = test-macros test-bisection \
	test-falsepos test-secant test-newton
EXTRA_DIST = test-macros test-bisection test-falsepos test-secant test-newton

TESTS_ENVIRONMENT = BUILDDIR=@full_top_builddir@/roots SRCDIR=@full_top_srcdir@ $(SHELL) -x

test_LDADD = libgslroots.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslroots.a

test_SOURCES = test.c test.h

libgslroots_a_SOURCES = utility.c bisection.c falsepos.c secant.c newton.c
