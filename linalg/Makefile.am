noinst_LTLIBRARIES = libgsllinalg.la 

check_PROGRAMS = test_la

pkginclude_HEADERS = gsl_linalg.h

noinst_HEADERS = tridiag.h

INCLUDES= -I$(top_builddir)

TESTS = test_la

test_la_LDADD = libgsllinalg.la ../blas/libgslblas.la ../blas/libgslblasnative.la ../sys/libgslsys.la ../err/libgsltest.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../permutation/libgslpermutation.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../utils/libutils.la

test_la_DEPENDENCIES = libgsllinalg.la
test_la_SOURCES = test_la.c

libgsllinalg_la_SOURCES = multiply.c tridiag.c tridiag.h lu.c hh.c qr.c qrpt.c svd.c householder.c

noinst_HEADERS =  givens.c matrix.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
