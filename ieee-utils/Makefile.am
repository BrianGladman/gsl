noinst_LTLIBRARIES = libgslieeeutils.la 

pkginclude_HEADERS = gsl_ieee_utils.h

libgslieeeutils_la_SOURCES =  print.c make_rep.c gsl_ieee_utils.h env.c fp.c read.c

noinst_HEADERS = fp-aix.c fp-darwin.c fp-hpux.c fp-hpux11.c fp-irix.c fp-gnum68k.c fp-gnuppc.c fp-solaris.c fp-gnusparc.c fp-sunos4.c fp-tru64.c fp-unknown.c fp-gnux86.c fp-freebsd.c fp-os2emx.c fp-netbsd.c fp-openbsd.c fp-gnuc99.c endian.c standardize.cun

INCLUDES= -I$(top_builddir)

TESTS = test
check_PROGRAMS = test
test_LDADD = libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la
test_SOURCES = test.c
