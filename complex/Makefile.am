noinst_LTLIBRARIES = libgslcomplex.la 

pkginclude_HEADERS = gsl_complex.h gsl_complex_math.h

INCLUDES = -I$(top_builddir)

libgslcomplex_la_SOURCES = math.c

TESTS = test
check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslcomplex.la ../err/libgslerr.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../err/libgsltest.la ../sys/libgslsys.la  ../utils/libutils.la
test_DEPENDENCIES = libgslcomplex.la
