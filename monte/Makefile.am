# -*-makefile-*-
## Libtool version
##

pkglib_LIBRARIES = libgslmonte.a 

libgslmonte_a_SOURCES = vegas.c vegas-print.c vegas_print.h miser.c plain.c \
			gsl_monte_vegas.h gsl_monte_miser.h gsl_monte_plain.h \
			gsl_monte.h utils.c utils.h

EXTRA_libgslmonte_a_SOURCES = init_source.c

check_PROGRAMS = vegas_test miser_test plain_test

pkginclude_HEADERS = gsl_monte_vegas.h gsl_monte_miser.h gsl_monte_plain.h \
		  gsl_monte.h

#noinst_HEADERS = vegas.h
noinst_HEADERS = utils.h init_source.c vegas_print.h

INCLUDES = -I$(srcdir)/../ -I$(srcdir)/../err -I$(srcdir)/../rng/ \
		-I$(srcdir)/../vector/

TESTS = vegas_test miser_test plain_test

TESTS_ENVIRONMENT = #BUILDDIR=@full_top_builddir@/roots SRCDIR=@full_top_srcdir@ $(SHELL) -x

vegas_test_SOURCES = vegas_test.c gsl_monte_vegas.h gsl_monte.h

vegas_test_LDADD = libgslmonte.a ../rng/libgslrng.a ../vector/libgslvector.a \
		../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

miser_test_SOURCES = miser_test.c gsl_monte_miser.h gsl_monte.h

miser_test_LDADD = libgslmonte.a ../rng/libgslrng.a ../vector/libgslvector.a \
		../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

plain_test_SOURCES = plain_test.c gsl_monte_plain.h gsl_monte.h

plain_test_LDADD = libgslmonte.a ../rng/libgslrng.a ../vector/libgslvector.a \
		../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
