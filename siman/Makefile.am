## Process this file with automake to produce Makefile.in

check_PROGRAMS = test siman_test siman_tsp
pkglib_LIBRARIES = libgslsiman.a

TESTS = test
EXTRA_DIST = siman_test_driver.sh

test_SOURCES = test.c
test_LDADD = libgslsiman.a ../rng/libgslrng.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslsiman.a

siman_test_SOURCES = siman_test.c
siman_test_LDADD = libgslsiman.a ../rng/libgslrng.a ../err/libgslerr.a ../utils/libutils.a
siman_test_DEPENDENCIES = libgslsiman.a

siman_tsp_SOURCES = siman_tsp.c
siman_tsp_LDADD = libgslsiman.a ../rng/libgslrng.a ../err/libgslerr.a ../utils/libutils.a
siman_tsp_DEPENDENCIES = libgslsiman.a

CLEANFILES = siman_test.out1 siman_test.out2

libgslsiman_a_SOURCES = siman.c

pkginclude_HEADERS = gsl_siman.h

INCLUDES = -I$(srcdir) -I$(srcdir)/../rng -I$(srcdir)/../err
check_INCLUDES = -I. -I../rng -I../err
