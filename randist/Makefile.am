B:1,$
## $Id$
## Assumes that ../random/gsl_ran.h and ../random/libgslran.a exist
## ...and are built!

bin_PROGRAMS=exponential gamma gauss poisson 
noinst_PROGRAMS=bench-gamma bench-gauss bench-poisson  \
   test-gauss test-shuffle

include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../random
LDADD= ../random/libgslrandom.a

gauss_SOURCES= rungauss.c gauss.c
test_gauss_SOURCES= test-gauss.c gauss.c
bench_gauss_SOURCES= bench-gauss.c gauss.c

exponential_SOURCES= runexponential.c exponential.c

poisson_SOURCES= runpoisson.c poisson.c 
poisson_LDADD= exponential.o $(LDADD)
bench_poisson_SOURCES= bench-poisson.c poisson.c
bench_poisson_LDADD= $(poisson_LDADD)

gamma_SOURCES= rungamma.c gamma.c
bench_gamma_SOURCES= bench-gamma.c gamma.c

test_shuffle_SOURCES= test-shuffle.c shuffle.c




