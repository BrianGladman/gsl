lib_LIBRARIES = libgslrng.a 

include_HEADERS = gsl_rng.h

INCLUDES = -I$(srcdir)/..  -I$(srcdir)/../err

libgslrng_a_SOURCES = gsl_rng.h rng.c default.c \
bsdrand.c cmrg.c minstd.c mrg.c mt.c r250.c ran0.c \
ran1.c ran2.c ran3.c rand.c rand48.c randu.c ranf.c ranlux.c ranmar.c \
taus.c tds.c tmp.c tt.c uni.c uni32.c vax.c zuf.c

test_SOURCES = test.c
test_LDADD = libgslrng.a ../err/libgslerr.a ../err/libgsltest.a
TESTS = test
check_PROGRAMS = test

#benchmark_PROGRAMS = benchmark
#benchmark_SOURCES = benchmark.c 
#benchmark_LDADD = libgslrng.a ../err/libgslerr.a