lib_LIBRARIES = libgslrandist.a
libgslrandist_a_SOURCES = \
        beta.c bernoulli.c binomial.c cauchy.c chisq.c erlang.c \
        exponential.c fdist.c flat.c \
	gamma.c gauss.c geometric.c levy.c logistic.c lognormal.c \
	nbinomial.c \
	pareto.c poisson.c shuffle.c sphere.c tdist.c twosidedexp.c weibull.c

include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../rng -I$(srcdir)/../specfunc -I$(srcdir)/../err


#bin_PROGRAMS = gsl-dist

#gsl_dist_SOURCES = gsl-dist.c
#gsl_dist_LDADD = libgslrandist.a ../rng/libgslrng.a ../specfunc/libgslspecfunc.a ../err/libgslerr.a

test_SOURCES = test.c
test_LDADD = libgslrandist.a ../rng/libgslrng.a ../specfunc/libgslspecfunc.a ../err/libgslerr.a ../err/libgsltest.a
TESTS = test
check_PROGRAMS = test
