# -*-makefile-*-
## Libtool version
##

noinst_LTLIBRARIES = libgslmonte.la 

libgslmonte_la_SOURCES = vegas.c vegas-print.c vegas_print.h miser.c plain.c \
			gsl_monte_vegas.h gsl_monte_miser.h gsl_monte_plain.h \
			gsl_monte.h utils.c utils.h

EXTRA_libgslmonte_la_SOURCES = init_source.c

check_PROGRAMS = vegas_test miser_test plain_test

pkginclude_HEADERS = gsl_monte_vegas.h gsl_monte_miser.h gsl_monte_plain.h \
		  gsl_monte.h

#noinst_HEADERS = vegas.h
noinst_HEADERS = utils.h init_source.c vegas_print.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)


TESTS = vegas_test miser_test plain_test

TESTS_ENVIRONMENT = #$(SHELL) -x

vegas_test_SOURCES = vegas_test.c gsl_monte_vegas.h gsl_monte.h

vegas_test_LDADD = libgslmonte.la ../rng/libgslrng.la ../vector/libgslvector.la \
	../block/libgslblock.la	../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la

miser_test_SOURCES = miser_test.c gsl_monte_miser.h gsl_monte.h

miser_test_LDADD = libgslmonte.la ../rng/libgslrng.la ../vector/libgslvector.la \
		../block/libgslblock.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la

plain_test_SOURCES = plain_test.c gsl_monte_plain.h gsl_monte.h

plain_test_LDADD = libgslmonte.la ../rng/libgslrng.la ../vector/libgslvector.la \
		../block/libgslblock.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.la


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
