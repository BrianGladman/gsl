noinst_LTLIBRARIES = libgslpoly.la 

pkginclude_HEADERS = gsl_poly.h

INCLUDES= -I$(top_builddir)

libgslpoly_la_SOURCES = solve_quadratic.c solve_cubic.c \
	zsolve_quadratic.c zsolve_cubic.c zsolve.c zsolve_init.c

noinst_HEADERS = balance.c companion.c norm.c qr.c

TESTS = test

check_PROGRAMS = test

test_SOURCES = test.c test_poly.c test.h
test_LDADD = libgslpoly.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la
test_DEPENDENCIES = libgslpoly.la
