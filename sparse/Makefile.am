noinst_LTLIBRARIES = libgslsparse.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_sparse.h

libgslsparse_la_SOURCES = spcompress.c spcopy.c spdgemm.c spdgemv.c spgetset.c spmatrix.c spoper.c spprop.c spswap.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslsparse.la ../test/libgsltest.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la  ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../sys/libgslsys.la ../err/libgslerr.la ../utils/libutils.la ../rng/libgslrng.la

test_SOURCES = test.c
