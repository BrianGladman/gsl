lib_LTLIBRARIES = libgslieeeutils.la 

include_HEADERS = gsl_ieee_utils.h

libgslieeeutils_la_SOURCES = endian.c print-ieee.c standardize.c \
	make_rep.c gsl_ieee_utils.h ieee_utils.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

TESTS = test
check_PROGRAMS = test
test_LDADD = libgslieeeutils.la ../err/libgsltest.a
test_SOURCES = test.c
