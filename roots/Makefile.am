# -*-makefile-*-

noinst_LIBRARIES = libgslroots.a 

check_PROGRAMS = test

include_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

TESTS = test

test_SOURCES = test.c test_macros.c test_funcs.c test_roots.c test_poly.c test.h
test_LDADD = libgslroots.a ../err/libgslerr.a ../err/libgsltest.a
test_DEPENDENCIES = libgslroots.a

libgslroots_a_SOURCES = bisection.c brent.c falsepos.c secant.c newton.c \
  solve_quadratic.c solve_cubic.c zsolve_quadratic.c zsolve_cubic.c
