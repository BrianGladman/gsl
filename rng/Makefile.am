lib_LIBRARIES = libgslrng.a 

include_HEADERS = gsl_rng.h

INCLUDES = -I$(srcdir)/..  -I$(srcdir)/../err

libgslrng_a_SOURCES = gsl_rng.h rng.c \
                      cmrg.c mrg.c rand.c taus.c uni.c uni32.c zuf.c \
		      minstd.c mt.c tt.c r250.c ran0.c ran1.c ran2.c ran3.c \
                      bsdrand.c randu.c vax.c ranlux.c \
                      default.c

test_SOURCES = test.c
test_LDADD = libgslrng.a ../err/libgslerr.a ../err/libgsltest.a
TESTS = test
check_PROGRAMS = test benchmark

benchmark_SOURCES = benchmark.c 
benchmark_LDADD = libgslrng.a ../err/libgslerr.a