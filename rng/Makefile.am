noinst_LIBRARIES = libgslrng.a 

include_HEADERS = gsl_rng.h

INCLUDES = -I$(srcdir)/..  -I$(srcdir)/../err

libgslrng_a_SOURCES = gsl_rng.h rng.c default.c \
cmrg.c slatec.c minstd.c mrg.c mt.c r250.c ran0.c \
ran1.c ran2.c ran3.c rand.c random.c rand48.c randu.c \
ranf.c ranlux.c ranmar.c \
taus.c transputer.c tt.c uni.c uni32.c vax.c zuf.c

test_SOURCES = test.c
test_LDADD = libgslrng.a ../err/libgslerr.a ../err/libgsltest.a
TESTS = test
check_PROGRAMS = test benchmark

benchmark_SOURCES = benchmark.c 
benchmark_LDADD = libgslrng.a ../err/libgslerr.a