lib_LIBRARIES = libgslspecfunc.a

# check_PROGRAMS = test_blas

include_HEADERS = \
	gsl_blas.h		\
	gsl_blas_types.h	\
	gsl_blas_raw.h

# noinst_HEADERS = 

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../vector -I$(srcdir)/../matrix

# TESTS = test_blas

#test_blas_LDADD = libgslblas.a ../err/libgslerr.a ../err/libgsltest.a
#test_blas_DEPENDENCIES = libgslblas.a
#test_blas_SOURCES = test_blas.c

libgslspecfunc_a_SOURCES = 	\
	blas.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
