# -*-makefile-*-
## Libtool version
##

lib_LIBRARIES = libgslmonte.a 

libgslmonte_a_SOURCES = vegas.c vegas-print.c miser.c plain.c gsl_vegas.h gsl_miser.h gsl_monte_plain.h gsl_monte.h

check_PROGRAMS = vegas_test miser_test plain_test

include_HEADERS = gsl_vegas.h gsl_miser.h gsl_monte_plain.h gsl_monte.h

#noinst_HEADERS = vegas.h

INCLUDES = -I$(srcdir)/../ -I$(srcdir)/../err -I$(srcdir)/../rng/ \
		-I$(srcdir)/../vector/

TESTS = vegas_test miser_test plain_test

TESTS_ENVIRONMENT = #BUILDDIR=@full_top_builddir@/roots SRCDIR=@full_top_srcdir@ $(SHELL) -x

vegas_test_SOURCES = vegas_test.c gsl_vegas.h gsl_monte.h

vegas_test_LDADD = libgslmonte.a ../rng/libgslrng.a ../vector/libgslvector.a \
		../err/libgslerr.a ../err/libgsltest.a

vegas_test_LDFLAGS = -static

miser_test_SOURCES = miser_test.c gsl_miser.h gsl_monte.h

miser_test_LDADD = libgslmonte.a ../rng/libgslrng.a ../vector/libgslvector.a \
		../err/libgslerr.a ../err/libgsltest.a 

miser_test_LDFLAGS = -static

plain_test_SOURCES = plain_test.c gsl_plain.h gsl_monte.h

plain_test_LDADD = libgslmonte.a ../rng/libgslrng.a ../vector/libgslvector.a \
		../err/libgslerr.a ../err/libgsltest.a 

plain_test_LDFLAGS = -static


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
