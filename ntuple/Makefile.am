noinst_LTLIBRARIES = libgslntuple.la

pkginclude_HEADERS = gsl_ntuple.h

INCLUDES= -I$(top_builddir)

libgslntuple_la_SOURCES = ntuple.c

test_SOURCES = test.c
test_LDADD = libgslntuple.la ../histogram/libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

TESTS = test
check_PROGRAMS = test

CLEANFILES = test.dat