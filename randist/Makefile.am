## Assumes that ../random/gsl_ran.h and ../random/libgslran.a exist
## ...land are built!

bin_PROGRAMS=gsl-exponential gsl-gamma gsl-gauss gsl-poisson \
	gsl-lorentz gsl-laplace
noinst_PROGRAMS= test-gauss test-shuffle

lib_LIBRARIES = libgslrandist.a #libutils.a
libgslrandist_a_SOURCES = gauss.c poisson.c gamma.c exponential.c lorentz.c \
	laplace.c

include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../random
LDADD= ../random/libgslrandom.a ../utils/libutils.a

gsl_gauss_SOURCES= gsl-gauss.c gauss.c
gsl_exponential_SOURCES= gsl-exponential.c exponential.c
gsl_lorentz_SOURCES = gsl-lorentz.c lorentz.c
gsl_laplace_SOURCES = gsl-laplace.c laplace.c
gsl_poisson_SOURCES= gsl-poisson.c poisson.c 
gsl_poisson_LDADD= exponential.o $(LDADD)
gsl_gamma_SOURCES= gsl-gamma.c gamma.c

test_gauss_SOURCES= test-gauss.c gauss.c
test_shuffle_SOURCES= test-shuffle.c shuffle.c
