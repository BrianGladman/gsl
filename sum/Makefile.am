pkglib_LIBRARIES = libgslsum.a 

check_PROGRAMS = test

pkginclude_HEADERS = \
	gsl_sum.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err 

TESTS = test

test_LDADD = libgslsum.a ../err/libgslerr.a ../err/libgsltest.a
test_DEPENDENCIES = libgslsum.a
test_SOURCES = test.c

libgslsum_a_SOURCES = 	\
	levin_u.c levin_uerr.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
