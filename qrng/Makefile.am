noinst_LTLIBRARIES = libgslqrng.la 

pkginclude_HEADERS = gsl_qrng.h


INCLUDES= -I$(top_builddir)

libgslqrng_la_SOURCES = \
	gsl_qrng.h	\
	qrng.c 		\
	niederreiter-2.c


test_SOURCES = test.c
test_LDADD = libgslqrng.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la

TESTS = test
check_PROGRAMS = test
