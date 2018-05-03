noinst_LTLIBRARIES = libgslfilter.la

pkginclude_HEADERS = gsl_filter.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslfilter_la_SOURCES = gaussian.c impulse.c median.c rmedian.c

noinst_HEADERS = test_gaussian.c test_rmedian.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslfilter.la ../movstat/libgslmovstat.la ../statistics/libgslstatistics.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../randist/libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../complex/libgslcomplex.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la ../vector/libgslvector.la ../block/libgslblock.la
