lib_LIBRARIES = libgslrandist.a
libgslrandist_a_SOURCES = cauchy.c chisq.c erlang.c exponential.c flat.c \
	gamma.c gauss.c geometric.c laplace.c logistic.c lognormal.c \
	pareto.c poisson.c shuffle.c weibull.c

include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../rng
LDADD= ../rng/libgslrng.a

