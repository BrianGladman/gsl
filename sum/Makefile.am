noinst_LTLIBRARIES = libgslsum.la 

check_PROGRAMS = test

pkginclude_HEADERS = \
	gsl_sum.h

INCLUDES= -I$(top_builddir)

TESTS = test

test_LDADD = libgslsum.la ../ieee-utils/libgslieeeutils.la ../sys/libgslsys.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la
test_DEPENDENCIES = libgslsum.la
test_SOURCES = test.c

libgslsum_la_SOURCES = 	\
	levin_u.c levin_uerr.c

