noinst_LIBRARIES = libgsllinalg.a 

check_PROGRAMS = test_la

include_HEADERS = \
	gsl_linalg.h

noinst_HEADERS = 


INCLUDES = \
	-I$(srcdir)/..		\
	-I$(srcdir)/../err	\
	-I$(srcdir)/../matrix	\
	-I$(srcdir)/../vector

TESTS = test_la

test_la_LDADD = libgsllinalg.a 	\
	../sys/libgslsys.a	\
	../err/libgslerr.a	\
	../err/libgsltest.a

test_la_DEPENDENCIES = libgsllinalg.a
test_la_SOURCES = test_la.c

libgsllinalg_a_SOURCES = 	\
	linalg_simple.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
