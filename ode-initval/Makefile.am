pkglib_LIBRARIES = libgslodeiv.a 

check_PROGRAMS = test_odeiv

pkginclude_HEADERS = \
	gsl_odeiv.h

noinst_HEADERS = odeiv_util.h


INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../ieee-utils -I$(srcdir)/../err -I$(srcdir)/../linalg -I$(srcdir)/../permutation -I$(srcdir)/../matrix -I$(srcdir)/../vector -I$(srcdir)/../block 

TESTS = test_odeiv

test_odeiv_LDADD = libgslodeiv.a \
	../ieee-utils/libgslieeeutils.a \
	../sys/libgslsys.a \
	../err/libgslerr.a \
	../err/libgsltest.a \
	../linalg/libgsllinalg.a \
	../matrix/libgslmatrix.a \
	../vector/libgslvector.a \
	../block/libgslblock.a \
        ../utils/libutils.a
test_odeiv_DEPENDENCIES = libgslodeiv.a
test_odeiv_SOURCES = test_odeiv.c

libgslodeiv_a_SOURCES = \
	control.c	\
	evolve.c	\
	monitor.c	\
	bsimp.c		\
	gear1.c		\
	gear2.c		\
	odeiv.c		\
	odeiv_util.c	\
	rk2.c		\
	rk4.c		\
	rk2imp.c	\
	rk4imp.c	\
	rk8pd.c		\
	rkck.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
