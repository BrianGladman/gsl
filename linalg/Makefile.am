noinst_LTLIBRARIES = libgsllinalg.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_linalg.h

noinst_HEADERS = tridiag.h

INCLUDES= -I$(top_builddir)

TESTS = test

test_LDADD = libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../sys/libgslsys.la ../test/libgsltest.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../permutation/libgslpermutation.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../utils/libutils.la

test_DEPENDENCIES = libgsllinalg.la
test_SOURCES = test.c

libgsllinalg_la_SOURCES = multiply.c tridiag.c tridiag.h lu.c hh.c qr.c qrpt.c svd.c householder.c cholesky.c

noinst_HEADERS =  givens.c

