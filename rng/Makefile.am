lib_LIBRARIES = libgslrng.a 

include_HEADERS = gsl_rng.h

INCLUDES = -I$(srcdir)/..  -I$(srcdir)/../err

libgslrng_a_SOURCES = \
	gsl_rng.h rng.c blah.c cmrg.c

test_SOURCES = test.c
test_LDADD = libgslrng.a ../err/libgslerr.a ../err/libgsltest.a
TESTS = test
check_PROGRAMS = test 
