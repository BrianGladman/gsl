## $Id$
## Assumes that ../ranunif/gsl_ran.h and ../ranunif/libgslran.a exist
## ...and are built!

bin_PROGRAMS= gauss test-gauss bench-gauss poisson bench-poisson test-shuffle gamma bench-gamma

include_HEADERS= gsl_randist.h

INCLUDES= -I../ranunif
LDADD= ../ranunif/libgslrantaus.a

gauss_SOURCES= rungauss.c gauss.c
test_gauss_SOURCES= test-gauss.c gauss.c
bench_gauss_SOURCES= bench-gauss.c gauss.c

poisson_SOURCES= runpoisson.c poisson.c 
poisson_LDADD= exponential.o $(LDADD)
bench_poisson_SOURCES= bench-poisson.c poisson.c
bench_poisson_LDADD= $(poisson_LDADD)

gamma_SOURCES= rungamma.c gamma.c
bench_gamma_SOURCES= bench-gamma.c gamma.c

test_shuffle_SOURCES= test-shuffle.c shuffle.c




