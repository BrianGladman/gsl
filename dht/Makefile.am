noinst_LTLIBRARIES = libgsldht.la 

check_PROGRAMS = test_dht

pkginclude_HEADERS = \
	gsl_dht.h

INCLUDES= -I$(top_builddir)

TESTS = test_dht

test_dht_LDADD = libgsldht.la		\
		../specfunc/libgslspecfunc.la	\
		../sys/libgslsys.la 	\
		../err/libgslerr.la	\
		../err/libgsltest.la     \
                ../utils/libutils.la

test_dht_DEPENDENCIES = libgsldht.la
test_dht_SOURCES = test_dht.c

libgsldht_la_SOURCES = 	\
		dht_transform.c


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
