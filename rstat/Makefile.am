noinst_LTLIBRARIES = libgslrstat.la

pkginclude_HEADERS = gsl_rstat.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslrstat_la_SOURCES = rstat.c rquantile.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslrstat.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la ../vector/libgslvector.la ../rng/gslrng.la ../randist/gslrandist.la


