lib_LIBRARIES = libgslieeeutils.a 

include_HEADERS = gsl_ieee_utils.h

libgslieeeutils_a_SOURCES = endian.c print-ieee.c standardize.c \
	make_rep.c gsl_ieee_utils.h ieee_utils.h \
	fp-env.c fp-mode.c 

noinst_HEADERS = fp-hpux.c fp-linux.c fp-unknown.c

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

TESTS = test
check_PROGRAMS = test
test_LDADD = libgslieeeutils.a ../err/libgsltest.a
test_SOURCES = test.c
