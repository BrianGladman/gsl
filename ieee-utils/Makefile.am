lib_LIBRARIES = libgslieeeutils.a 

include_HEADERS = gsl_ieee_utils.h

libgslieeeutils_a_SOURCES = endian.c print-ieee.c standardize.c \
	gsl_ieee_utils.h ieee_utils.h

check_PROGRAMS = test
test_LDADD = libgslieeeutils.a
test_SOURCES = test.c
