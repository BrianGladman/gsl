# -*-makefile-*-

lib_LTLIBRARIES = libgslroots.la 

check_PROGRAMS = test

include_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES = -I$(srcdir)/../err

TESTS = test-macros test-bisection \
	test-falsepos test-secant test-newton
EXTRA_DIST = test-macros test-bisection test-falsepos test-secant test-newton

TESTS_ENVIRONMENT = BUILDDIR=@full_top_builddir@/roots SRCDIR=@full_top_srcdir@ $(SHELL) -x

test_LDADD = libgslroots.la ../err/libgslerr.la ../err/libgsltest.la

test_DEPENDENCIES = libgslroots.la

test_SOURCES = test.c test.h

libgslroots_la_SOURCES = utility.c bisection.c falsepos.c secant.c newton.c
