lib_LIBRARIES = libgslrandist.a
libgslrandist_a_SOURCES = gauss.c poisson.c gamma.c exponential.c laplace.c \
			  cauchy.c flat.c

include_HEADERS= gsl_randist.h

INCLUDES= -I$(srcdir)/../rng
LDADD= ../rng/libgslrng.a

