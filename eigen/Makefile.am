noinst_LTLIBRARIES = libgsleigen.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_eigen.h
libgsleigen_la_SOURCES = eigen_sort.c jacobi.c
INCLUDES= -I$(top_builddir)

TESTS = test

test_LDADD = libgsleigen.la ../sys/libgslsys.la ../test/libgsltest.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../blas/libgslblasnative.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../err/libgslerr.la ../utils/libutils.la

test_DEPENDENCIES = libgsleigen.la
test_SOURCES = test.c


