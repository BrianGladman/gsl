pkglib_LIBRARIES = libgsllinalg.a 

check_PROGRAMS = test_la

pkginclude_HEADERS = \
	gsl_linalg.h

noinst_HEADERS = tridiag.h


INCLUDES = \
	-I$(srcdir)/..		\
	-I$(srcdir)/../err	\
	-I$(srcdir)/../matrix	\
	-I$(srcdir)/../vector

TESTS = test_la

test_la_LDADD = libgsllinalg.a 		\
	../sys/libgslsys.a		\
	../err/libgsltest.a		\
	../matrix/libgslmatrix.a	\
	../vector/libgslvector.a	\
	../err/libgslerr.a

test_la_DEPENDENCIES = libgsllinalg.a
test_la_SOURCES = test_la.c

libgsllinalg_a_SOURCES = 	\
	eigen_sort.c		\
	linalg_simple.c		\
	multiply.c		\
	tridiag.c


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
