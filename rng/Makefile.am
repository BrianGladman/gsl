lib_LIBRARIES = libgslrng.a 

include_HEADERS = gsl_rng.h

INCLUDES = -I$(srcdir)/..  -I$(srcdir)/../err

libgslrng_a_SOURCES = gsl_rng.h rng.c \
                      cmrg.c mrg.c rand.c taus.c uni.c uni32.c zuf.c \
		      minstd.c \
                      bad_rand.c bad_randu.c vax.c ranlux.c \
                      default.c

test_SOURCES = test.c
test_LDADD = libgslrng.a ../err/libgslerr.a ../err/libgsltest.a
TESTS = test
check_PROGRAMS = test 
