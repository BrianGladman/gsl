lib_LIBRARIES = libgslrandist.a
libgslrandist_a_SOURCES = \
        beta.c binomial.c cauchy.c chisq.c erlang.c \
        exponential.c fdist.c flat.c \
	gamma.c gauss.c geometric.c laplace.c logistic.c lognormal.c \
	pareto.c poisson.c shuffle.c sphere.c tdist.c weibull.c

include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../rng
LDADD= ../rng/libgslrng.a

