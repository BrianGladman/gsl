noinst_LTLIBRARIES = libgslodeiv.la 

pkginclude_HEADERS = gsl_odeiv.h

noinst_HEADERS = odeiv_util.h

INCLUDES= -I$(top_builddir)

libgslodeiv_la_SOURCES = control.c evolve.c monitor.c bsimp.c gear1.c gear2.c odeiv.c odeiv_util.c rk2.c rk4.c rk2imp.c rk4imp.c rk8pd.c rkck.c

check_PROGRAMS = test_odeiv

TESTS = test_odeiv

test_odeiv_LDADD = libgslodeiv.la ../ieee-utils/libgslieeeutils.la ../sys/libgslsys.la ../err/libgslerr.la ../err/libgsltest.la ../linalg/libgsllinalg.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../permutation/libgslpermutation.la ../block/libgslblock.la ../utils/libutils.la
test_odeiv_DEPENDENCIES = libgslodeiv.la
test_odeiv_SOURCES = test_odeiv.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
