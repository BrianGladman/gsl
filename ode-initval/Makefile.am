pkglib_LIBRARIES = libgslodeiv.a 

check_PROGRAMS = test_odeiv

pkginclude_HEADERS = \
	gsl_odeiv.h

noinst_HEADERS = 


INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../ieee-utils -I$(srcdir)/../err 

TESTS = test_odeiv

test_odeiv_LDADD = libgslodeiv.a \
	../ieee-utils/libgslieeeutils.a \
	../sys/libgslsys.a \
	../err/libgslerr.a \
	../err/libgsltest.a
test_odeiv_DEPENDENCIES = libgslodeiv.a
test_odeiv_SOURCES = test_odeiv.c

libgslodeiv_a_SOURCES = \
	gear1.c		\
	gear2.c		\
	odeiv.c		\
	odeiv_util.c	\
	rk2.c		\
	rk4.c		\
	rk4imp.c	\
	rk8pd.c		\
	rkck.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
