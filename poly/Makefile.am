pkglib_LIBRARIES = libgslpoly.a 

pkginclude_HEADERS = gsl_poly.h

INCLUDES= -I$(top_builddir)

libgslpoly_a_SOURCES = solve_quadratic.c solve_cubic.c \
	zsolve_quadratic.c zsolve_cubic.c zsolve.c zsolve_init.c

noinst_HEADERS = balance.c companion.c norm.c qr.c

TESTS = test

check_PROGRAMS = test

test_SOURCES = test.c test_poly.c test.h
test_LDADD = libgslpoly.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslpoly.a
