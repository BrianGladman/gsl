pkglib_LIBRARIES = libgsldht.a 

check_PROGRAMS = test_dht

pkginclude_HEADERS = \
	gsl_dht.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err  -I$(srcdir)/../specfunc

TESTS = test_dht

test_dht_LDADD = libgsldht.a		\
		../specfunc/libgslspecfunc.a	\
		../sys/libgslsys.a 	\
		../err/libgslerr.a	\
		../err/libgsltest.a     \
                ../utils/libutils.a

test_dht_DEPENDENCIES = libgsldht.a
test_dht_SOURCES = test_dht.c

libgsldht_a_SOURCES = 	\
		dht_transform.c


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
