pkglib_LIBRARIES = libgslieeeutils.a 

pkginclude_HEADERS = gsl_ieee_utils.h

libgslieeeutils_a_SOURCES = endian.c print.c standardize.c \
	make_rep.c gsl_ieee_utils.h ieee_utils.h \
	env.c fp.c read.c

noinst_HEADERS = fp-hpux.c fp-linux.c fp-sunos4.c fp-solaris.c fp-unknown.c

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

TESTS = test
check_PROGRAMS = test
test_LDADD = libgslieeeutils.a ../err/libgsltest.a
test_SOURCES = test.c
