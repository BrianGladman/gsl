noinst_LIBRARIES = libgslrandist.a
libgslrandist_a_SOURCES = \
 bernoulli.c beta.c binomial.c cauchy.c chisq.c erlang.c exponential.c \
 exppow.c fdist.c flat.c gamma.c gauss.c geometric.c \
 gumbel.c hyperg.c laplace.c levy.c logarithmic.c logistic.c lognormal.c \
 nbinomial.c pareto.c pascal.c poisson.c rayleigh.c \
 shuffle.c sphere.c tdist.c weibull.c 


include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../rng -I$(srcdir)/../specfunc -I$(srcdir)/../err -I$(top_srcdir)


bin_PROGRAMS = gsl-randist

gsl_randist_SOURCES = gsl-randist.c
gsl_randist_LDADD = libgslrandist.a ../rng/libgslrng.a ../specfunc/libgslspecfunc.a ../err/libgslerr.a ../utils/libutils.a ../sys/libgslsys.a

test_SOURCES = test.c
test_LDADD = libgslrandist.a ../rng/libgslrng.a ../specfunc/libgslspecfunc.a ../err/libgslerr.a ../err/libgsltest.a  ../utils/libutils.a ../sys/libgslsys.a
TESTS = test
check_PROGRAMS = test
