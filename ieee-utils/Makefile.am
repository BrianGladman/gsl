pkglib_LIBRARIES = libgslieeeutils.a 

pkginclude_HEADERS = gsl_ieee_utils.h

libgslieeeutils_a_SOURCES =  print.c \
	make_rep.c gsl_ieee_utils.h \
	env.c fp.c read.c

noinst_HEADERS = fp-hpux.c fp-linux.c fp-sunos4.c fp-solaris.c fp-tru64.c \
        fp-sparclinux.c fp-unknown.c endian.c standardize.c

INCLUDES= -I$(top_builddir)

TESTS = test
check_PROGRAMS = test
test_LDADD = libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_SOURCES = test.c
